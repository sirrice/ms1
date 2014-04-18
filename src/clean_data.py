import csv
import pdb
import numpy as np


def rewrite_full_data(infilename='clean10_4.txt', outfilename='out'):
  """
  Strips input of everything except M/Z and intensity values.
  Removes all readings where intensity is below 10

  Input format: 
    file of lists of numbers where odd lines are the m/z values and
    even lines are the corresponding intensities

    Each line has the format:

      binary: [NUMBER] <LIST OF NUMBERS>

  Output format:
    Similar to input format where each line has format:

      <LIST OF NUMBERS>
    
  """
  data = []
  with file(infilename) as f:
    scan_no = 0
    for l in f:
      arr = l.strip().split(" ")
      arr = arr[2:]
      xs = map(float, arr)

      l = f.next()
      arr = l.strip().split(" ")
      arr = arr[2:]
      ys = map(float, arr)

      xys = zip(xs, ys)
      subset = np.array(xys)

      # filter out intensities below 10
      subset = subset[subset[:,1] > 10]
      data.append(subset)


      scan_no += 1
      if scan_no % 10 == 0:
        print scan_no


  print "writing to %s" % outfilename
  with file(outfilename, 'w') as f:
    for subset in data:
      xs = subset[:,0]
      ys = subset[:,1]
      f.write(' '.join(map(str, xs)))
      f.write('\n')
      f.write(' '.join(map(str, ys)))
      f.write('\n')



# deprecated
def process_shitty_data():
  data = []
  with file('cleaned_ms.txt') as f:
    dialect = csv.Sniffer().sniff(f.read(1024))
    f.seek(0)
    reader = csv.reader(f, dialect)
    header = reader.next()
    print header

    for l in reader:
      try:
        l = map(float, l)
        data.append(l)
      except:
        print l
        exit()


  data = np.asarray(data)
  vector = data[:,1] == 1
  data = data[ vector ]
  scans = sorted(set(data[:,0]))

  for scan_id in scans:
    print scan_id
    vector = data[:,0] == scan_id
    subset = data[vector]
    mzidx = 9
    intensityidx = -1
    subset = subset[:,(mzidx, intensityidx)]

    if scan_id > 2000:
      process_scan(subset)


if __name__ == '__main__':


  import sys
  if len(sys.argv) < 3:
    print "python clean_data.py <input file name e.g., clean10_4.txt> <output file name e.g., out>"
    exit()
  infilename = sys.argv[1]
  outfilename = sys.argv[2]
  rewrite_full_data(infilename=infilename, outfilename=outfilename)
