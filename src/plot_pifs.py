import pdb
import csv
import time
import numpy as np
from collections import deque
from matplotlib import pyplot as plt


def load_pifs(fname='pifs'):
  """
  Read the output from process_all_data() and plot the data as a scatter plot
  """
  with file(fname) as f:
    data = np.array([map(float, l.strip().split(' ')) for l in f])

    print "scan\tminoldpif\tminnewpif\tnum pifs"
    for scan in sorted(set(data[:,0])):
      subset = data[data[:,0] == scan]
      idxs = np.arange(len(subset))
      old_pifs = subset[:,3]
      new_pifs = subset[:,-1]
      idxs = sorted(idxs, key=lambda idx: new_pifs[idx], reverse=True)
      idxs = idxs[:10]
      plt.scatter(old_pifs, new_pifs, lw=0, s=.75, alpha=.2, color='grey')

      olds, news = old_pifs[idxs], new_pifs[idxs]
      if scan % 100 == 0:
        print "%d\t%.4f -> %.4f\t%d" % (scan, np.percentile(olds, 50), np.percentile(news, 50), len(old_pifs))
      plt.scatter(olds, news, lw=0, alpha=.2, color='red')
    plt.show()
    pdb.set_trace()



def plot_old_news(fname='pifs'):
  """
  Read the output from process_all_data() and plot scan_no vs pif plots to compare before
  and after
  """
  with file(fname) as f:
    data = np.array([map(float, l.strip().split(' ')) for l in f])
    fig = plt.figure()
    ax1 = fig.add_subplot(4,1,1)
    ax2 = fig.add_subplot(4,1,2)
    ax3 = fig.add_subplot(4,1,3)
    ax4 = fig.add_subplot(4,1,4)
    prevs = deque()
    cum_improvments = []
    

    scan_nos = sorted(set(data[:,0]))
    for scan in scan_nos:
      subset = data[data[:,0] == scan]
      idxs = np.arange(len(subset))
      old_pifs = subset[:,3]
      new_pifs = subset[:,-1]
      xs = [scan] * len(old_pifs)

      # by intensity
      intensity_idxs = sorted(idxs, key=lambda idx: subset[idx,2], reverse=True)
      intensity_idxs = intensity_idxs[:10]
      ax1.scatter([scan]*len(intensity_idxs), old_pifs[intensity_idxs], 
        lw=0, s=1., alpha=.2, color='red')


      oldpif_idxs = sorted(idxs, key=lambda idx: subset[idx,3], reverse=True)
      oldpif_idxs = oldpif_idxs[:10]
      #oldpif_idxs = filter(lambda idx: not any([subset[idx,0] in seen for seen in prevs]), oldpif_idxs)
      ax2.scatter([scan]*len(oldpif_idxs), old_pifs[oldpif_idxs], 
        lw=0, s=3., alpha=.2, color='red')
      prevs.append(set(subset[oldpif_idxs,0]))
      if len(prevs) > 30:
        prevs.popleft()


      new_idxs = sorted(idxs, key=lambda idx: new_pifs[idx], reverse=True)
      new_idxs = new_idxs[:10]
      olds, news = old_pifs[new_idxs], new_pifs[new_idxs]
      #ax3.scatter(xs, new_pifs, lw=0, s=1., alpha=.2, color='grey')
      ax3.scatter([scan]*len(news), news, lw=0, s=3., alpha=.2, color='red')
      cum_improvments.append(sum(news) - sum(old_pifs[oldpif_idxs]))

      if scan % 100 == 0:
        print "%d\t%.4f -> %.4f\t%d" % (scan, np.percentile(olds, 50), np.percentile(news, 50), len(old_pifs))

    cum = 0.
    ys = []
    for cur in cum_improvments:
      cum += cur
      ys.append(cum)

    ax4.plot(scan_nos, ys, color='red')



    plt.show()
    pdb.set_trace()




if __name__ == '__main__':

  import sys
  if len(sys.argv) < 2:
    print "plot pifs"
    print "python plot_pifs.py <pif file>"
    exit()

  fname = sys.argv[1]
  plot_old_news(fname)
  #load_pifs(fname)
