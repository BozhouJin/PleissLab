import os, gzip, subprocess, HTSeq, sys
from itertools import zip_longest
from collections import defaultdict
import subprocess

def EmptyCounting(sample, counts, home_directory):

	with gzip.open("%s/noUmi/%s_R1.fastq.gz" % (home_directory, sample),'rt') as infile:
		
		

		for line in infile:
			
			if line[0] == '@':
				counts[(sample, 'total')] += 1
				
				info = line.rstrip().split(' ')
				seq = next(infile).rstrip()
				extra = next(infile).rstrip()
				quality = next(infile).rstrip()
				
				length = len(seq)

				if length <= 33:
					counts[(sample, 'empty')] += 1

			else:
				pass

def main():
		#start empty counting

	print("Empty Counting...")
	sample = 'wt1'

	counts = defaultdict(float)
	home_directory = '/home/jbz/Desktop/PleissLab/Illumina/4.122022/'
	EmptyCounting(sample, counts, home_directory)

	with open('/home/jbz/Desktop/PleissLab/Illumina/4.122022/EmptyCounts_20220426.txt', 'wt') as outfile:
		outfile.write('Sample\tTotal\tEmpty\tPercentEmpty\n')
		
		outfile.write('%s\t%f\t%f\n' % (sample, counts[(sample, 'total')], counts[(sample, 'empty')]))

main()