import pdb
import csv
import time
import numpy as np
from collections import deque
from matplotlib import pyplot as plt

dists = [.5, 1./3, .25]
mass = 1.00335483778 
dists = [1.00335483778 * dist for dist in dists]
dists = np.array(dists)


def eq(v1, v2, thresh=2e-5):
  return abs(v1-v2) < thresh

def is_isotope(x, cand, isotope_block, thresh=2e-5):
  """
  Optimized
  """
  diff = abs(x-cand)
  ppm = thresh * cand 
  if ppm > .1:
    pdb.set_trace()
  return diff <= 1.1*mass and (
      abs(diff) == 0 or
      abs(diff-isotope_block*mass) < ppm or
      abs(diff-2*isotope_block*mass) < ppm or
      abs(diff-3*isotope_block*mass) < ppm 
   )


def find_isotope_clusters(x, data):
  """
  There can be many fucking clusters defined by different
  m/z block sizes...

  Return:
    list of clusters of isotopes
  """
  xs = data[:,0]
  idxs = (xs >= x) & (xs <= x+(4*.51))
  xs = xs[idxs]
  clusters = []
  for isotope_block in [.25, 1./3, .5]:
    cluster_idxs = find_isotope(x, xs, isotope_block)
    if len(cluster_idxs) > 1:
      clusters.append(data[idxs][cluster_idxs,:])
  return clusters


def find_isotope(x, xs, isotope_block):
  """
  Find isotopes of x within data

  Args:
    x: m/z value
    data: np array of [[m/z, intensity], ...]

  Return:
    subset of data that are isotopes
  """
  retidxs = []
  for idx, cx in enumerate(xs):
    if is_isotope(x, cx, isotope_block):
      retidxs.append(idx)
  return retidxs

def in_window(bound, data):
  minx, maxx = tuple(bound)
  idxs = (data[:,0] >= minx) & (data[:,0] <= maxx)
  return data[idxs]

def compute_best_pif(data, x, isotopes, window):
  """
  return center of window containing best PIF score

  Args:
    x: current center of window
    isotopes: set of isotope x values.  used to compute numerator of pif 
  """
  minx, maxx = x - window, x + window
  window_data = in_window([minx, maxx], data)
  best_pif, best_x = None, None
  xs = window_data[:,0]
  ys = window_data[:,1]

  # start and end index of sliding window
  sidx, eidx = 0, 0

  # initialization
  isotope_intensity = 0.
  if window_data[0][0] in isotopes:
    isotope_intensity = window_data[0][1]
  total_intensity = window_data[0][1]


  for x, y in window_data:
    # if x extends outside max window size,
    # keep discarding left side of window
    while (x - xs[sidx]) > window/2.:
      if xs[sidx] in isotopes:
        isotope_intensity -= ys[sidx]
      total_intensity -= ys[sidx]
      sidx += 1

    # ??? /wu
    while eidx+1 < len(xs) and (xs[eidx+1] - x) <= window/2.:
      eidx += 1
      if xs[eidx] in isotopes:
        isotope_intensity += ys[eidx]
      total_intensity += ys[eidx]

    if total_intensity == 0: continue

    pif = isotope_intensity / total_intensity
    if best_pif is None or pif > best_pif:
      best_pif, best_x = pif, x

  return best_x, best_pif


def compute_centered_pif(data, x, isotopes, window):
  """
  compute the pif value for a window centered around x
  """
  minx, maxx = x - window/2, x + window/2
  window_data = in_window([minx, maxx], data)
  isotopes = in_window([minx, maxx], isotopes)
  isotope_intensity = sum(isotopes[:,1])
  total_intensity = sum(window_data[:,1])
  pif = isotope_intensity / total_intensity
  return pif



def compute_pifs(data, valid_indices, window=1.1):
  """
  Compute the windows that 
  * contains an actual peak by looking for isotopes
  * have the maximum possible PIF given all the other m/z 
    readings and the isotopes

  Args
    data: (m/z, intensity) pairs
    valid_indicies: index vals into data of the peaks that should
      be considered as isotopes
    window: window width

  Return
    List of (m/z value, intensity, pif, best m/z value to center window, best pif)
  """
  ret = []
  seen = set()
  for valid_idx in valid_indices:
    x, y = data[valid_idx]

    if x not in seen: 
      isotope_clusters = find_isotope_clusters(x, data)
      for isotope_cluster in isotope_clusters:
        ixs = isotope_cluster[:,0]
        seen.update(ixs)
        max_x, max_y = tuple(max(isotope_cluster, key=lambda p: p[1]))
        naive_pif = compute_centered_pif(data, max_x, isotope_cluster, 2.)
        best_x, best_pif = compute_best_pif(data, x, isotope_cluster, 2.)

        if best_pif is not None:
          best = (x, y, naive_pif, best_x, best_pif)
          ret.append(best)

  return ret

  # naive exhaustive search
  for x, y in data:
    minx, maxx = x - window/2, x + window/2
    window_data = in_window([minx, maxx], data)
    isotopes = find_isotope(x, window_data)

    if len(isotopes) == 0: continue
    if x < max(isotopes[:,0]): continue

    best = compute_pif_windows(data, x, y, isotopes, window)
    if best: ret.append(best)

  return filter(bool, ret)


def compute_valid_indices(data, prev_scans):
  """
  compute the indices of m/z values that should be considered as 
  potential isotopes (not noise)

  Return
    numpy array of the index valuse into data
  """
  ys = data[:,1]
  thresh = np.percentile(ys, 50)
  
  truncate = lambda x: round(x, 3)
  keys = set()
  for idx, prev_scan in enumerate(prev_scans):
    cur_xs = prev_scan[:,0]
    cur_xs = map(truncate, cur_xs)
    if idx == 0:
      keys = set(cur_xs)
    else:
      keys.intersection_update(set(cur_xs))

  idxs = [] #np.zeros(len(data)).astype(bool)
  for idx, x in enumerate(data[:,0]):
    if truncate(x) in keys:
      idxs.append(idx)

  return np.array(idxs)

def process_scan(scan_no, data, prev_scans):
  """
  look for peaks in the scan data
  """
  if len(prev_scans) < 2: return []
  if scan_no < 600: return []
  idxs = compute_valid_indices(data, prev_scans)
  if len(idxs) == 0: return []
  peaks = compute_pifs(data, idxs)

  if False and scan_no == 700:
    subset = data[idxs,:]
    xs = subset[:,0]
    xs = subset[:,0]
    ys = subset[:,1]

    #subset = subset[(xs >= 1240) & (xs <= 1260)]
    #xs, ys = subset[:,0], subset[:,1]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    thresh = np.percentile(subset[:,1], 29)
    for x,y in subset:
      ax.plot([x,x], [0, min(thresh, y)], alpha=0.3, color='grey')

    for x,y, pif, best_x, best_pif in peaks:
      y = min(thresh, y)
      ax.plot([best_x,best_x], [0, y], color='red')
      ax.plot([best_x-.5*mass,best_x-.5*mass], [.8*thresh, thresh], color='green')
      ax.plot([best_x+.5*mass,best_x+.5*mass], [.8*thresh, thresh], color='green')
    ax.set_xlim(450, 458)
    plt.show()

    pdb.set_trace()
  return peaks



def process_all_data(infname='out', outfname='pifs', stopping=-1):
  """
  process the output of rewrite_full_data()

  input file has format:

      M/Z values
      Intensity values
      M/Z values
      Intensity values

  Where each line is 

      <LIST OF FLOATS>


  """
  start = time.time()
  with file(outfname, 'w') as outf:
    with file(infname) as inf:
      q = deque()
      scan_no = 0
      for l in inf:
        xs = map(float, l.strip().split(' '))
        ys = map(float, inf.next().strip().split(' '))
        xys = zip(xs, ys)
        xys.sort()
        xys = np.array(xys)
        scan_no += 1

        if stopping > -1 and scan_no > stopping:
          print "stopping at scan %d" % scan_no
          break

        pif_list = process_scan(scan_no, xys, q)

        if scan_no % 10 == 0:
          print scan_no, len(pif_list)

        for tup in pif_list:
          row = [scan_no]
          row.extend(tup)
          outf.write("%s\n" % ' '.join(map(str, row)))

        if len(q) >= 2:
          q.popleft()
        q.append(xys)

  end = time.time()
  print "cost:\t", end-start,'\tperscan:', float(end-start)/scan_no


if __name__ == '__main__':

  import sys
  if len(sys.argv) < 3:
    print "processes cleaned MS1 data file and writes list of pifs"
    print "python clean_data.py <input file name e.g., out> <output file name e.g., pifs> <num scans to compute>"
    exit()
  infilename = sys.argv[1]
  outfilename = sys.argv[2]
  stopping = -1
  if len(sys.argv) > 3:
    stopping = int(sys.argv[3])



  process_all_data(infilename, outfilename, stopping)