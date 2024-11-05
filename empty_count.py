import os, gzip, subprocess, HTSeq, sys
from itertools import zip_longest
from collections import defaultdict
import subprocess

def EmptyCounting(sample, counts, home_directory):

	with gzip.open("%s/trimmed/%s_r1.fastq.gz" % (home_directory, sample),'rt') as infile:
		
		currentline = next(infile,None)
		while currentline != None:
			assert currentline[0] == '@'
			counts[(sample, 'total')] += 1	
			info = currentline.rstrip().split(' ')
			seq = next(infile).rstrip()
			extra = next(infile).rstrip()
			quality = next(infile).rstrip()	
			length = len(seq)

			if length <= 40:
				counts[(sample, 'empty')] += 1
			currentline = next(infile,None)


def main():
		#start empty counting

	print("Empty Counting...")
	home_directory = sys.argv[1]
	samples = sys.argv[2:len(sys.argv)]
	counts = defaultdict(float)
	for sample in samples:
		print("Counting %s"%(sample))
		EmptyCounting(sample, counts, home_directory)

	with open('%s/EmptyCounts.txt'%home_directory, 'wt') as outfile:
		outfile.write('Sample\tTotal\tEmpty\tPercentEmpty\n')
		for sample in samples:
			outfile.write('%s\t%f\t%f\t%f\n' % (sample, counts[(sample, 'total')], counts[(sample, 'empty')],(counts[(sample, 'empty')]/counts[(sample, 'total')])*100))

main()