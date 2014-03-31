import pdb
import csv
import numpy as np
from matplotlib import pyplot as plt

def eq(v1, v2, thresh=1e-3):
  return abs(v1-v2) < thresh

def is_isotope(x, cand):
  dists = [-.5, -1./3, -.25, .5, 1./3, .25]
  dists = [1.008 * dist for dist in dists]
  if abs(x-cand) > .51: return False

  for dist in dists:
    if eq(x+dist, cand):
      return True
  return False

def find_isotope(x, data):
  """
  Find isotopes of x within data

  Args:
    x: m/z value
    data: np array of [[m/z, intensity], ...]

  Return:
    subset of data that are isotopes
  """
  isotopes = []
  for cx, cy in data:
    if is_isotope(x, cx):
      isotopes.append((cx, cy))
  return np.array(isotopes)

def in_window(bound, data):
  minx, maxx = tuple(bound)
  vector = (data[:,0] >= minx) & (data[:,0] <= maxx)
  return data[vector]


def best_window(x, data, isotopes, window=1., total_steps=100.):
  """
  return center of new window using naive exhaustive search
  """
  minx, maxx = x - window, x + window
  window_data = in_window([minx, maxx], data)
  best_pif, best_x = None, None

  # For each possible window, compute its pif
  step_size = (maxx - minx) / total_steps
  for step_no in xrange(int(total_steps)):
    cur_x = (minx + window/2) + step_no * step_size
    minx, maxx = cur_x - window/2., cur_x + window/2.
    cur_data = in_window([minx, maxx], data)

    isotope_intensity = 0
    for x,y in cur_data:
      if x in isotopes:
        isotope_intensity += y
    total_intensity = sum(cur_data[:,1])
    if total_intensity == 0: continue

    pif = isotope_intensity / total_intensity
    if best_pif is None or pif > best_pif:
      best_pif, best_x = pif, cur_x

  return best_x, best_pif

def compute_pifs(data, window=1.):
  """
  Check each M/Z value if
  * it has isotopes within window
  * is the highest intensity isotope 
  * looks for best pif window using best_window()

  Return
    List of (m/z value, intensity, pif, best m/z value to center window, best pif)
  """
  ret = []
  for x, y in data:
    minx, maxx = x - window/2, x + window/2
    window_data = in_window([minx, maxx], data)

    isotopes = find_isotope(x, window_data)
    if len(isotopes) == 0: continue
    if x < max(isotopes[:,0]): continue

    best_x, best_pif = best_window(x, data, set(isotopes[:,0]), total_steps=20)
    if best_pif < .2: continue

    isotope_intensity = sum(isotopes[:,1])
    total_intensity = sum(window_data[:,1])
    pif = isotope_intensity / total_intensity
    if False and x != best_x:
      print x, best_x, pif, best_pif, '\t', y
    ret.append((x,y,pif, best_x, best_pif))


  return ret


def filter_data(data):
  thresh = np.percentile(data[:,1], 50)
  return data[data[:,1] > thresh]

def process_scan(scan_no, subset):
  """
  look for peaks in the scan data
  """
  subset = filter_data(subset)
  peaks = compute_pifs(subset)
  return peaks

  if False and peaks:
    xs = subset[:,0]
    ys = subset[:,1]

    #subset = subset[(xs >= 1240) & (xs <= 1260)]
    #xs, ys = subset[:,0], subset[:,1]
    plt.clf()
    thresh = np.percentile(subset[:,1], 99)
    for x,y in subset:
      plt.plot([x,x], [0, min(thresh, y)], alpha=0.3, color='grey')

    for x,y, pif in peaks:
      plt.plot([x,x], [0, y], color='red')
    plt.show()

    pdb.set_trace()



def process_all_data(fname='out'):
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
  all_pifs = []
  with file(fname) as f:
    scan_no = 0
    for l in f:
      xs = map(float, l.strip().split(' '))
      ys = map(float, f.next().strip().split(' '))
      xys = zip(xs, ys)
      xys.sort()
      xys = np.array(xys)
      scan_no += 1


      if scan_no % 10 == 0:
        print scan_no

      pif_list = process_scan(scan_no, xys)
      for tup in pif_list:
        row = [scan_no]
        row.extend(tup)
        all_pifs.append(row)

  # Write out pifs
  with file('pifs', 'w') as f:
    for row in all_pifs:
      f.write("%s\n" % ' '.join(map(str, row)))



def load_pifs(fname='pifs'):
  with file('pifs') as f:
    data = np.array([map(float, l.strip().split(' ')) for l in f])

    for scan in sorted(set(data[:,0])):
      subset = data[data[:,0] == scan]
      idxs = np.arange(len(subset))
      old_pifs = subset[:,3]
      new_pifs = subset[:,-1]
      idxs = sorted(idxs, key=lambda idx: new_pifs[idx], reverse=True)
      idxs = idxs[:10]
      plt.scatter(old_pifs, new_pifs, lw=0, alpha=.3, color='grey')

      olds, news = old_pifs[idxs], new_pifs[idxs]
      print scan, '\t', min(news), len(old_pifs)
      plt.scatter(olds, news, lw=0, alpha=.3, color='red')
    plt.show()
    pdb.set_trace()