import os, gzip, subprocess, sys, HTSeq
from itertools import zip_longest
from collections import defaultdict
import subprocess

def main():

	directory = sys.argv[1]
	samples = sys.argv[2:len(sys.argv)]
	for sample in samples:
		print("Empty filtering for %s"%(sample))
		EmptyFilter('%s/trimmed/'%(directory),'%s/noUmiNoEmpty/'%(directory),sample,33)

# This script is the same as remove_umi_filter_empty.py except for that it does not remove umi. This
# is specifically designed for my SPE-seq protocol since it does not contain the UMI in the primer.
# It removes empty reads, which as explained in remove_umi_filter_empty.py, has at most 26bp primer
# sequence, plus ~7bp N9 sequence (assuming the min length for N9 to anneal is 2bp), plus the 7bp
# UMI. 
# @param n: the maximum length of empty to be filtered out. Default is 33 for all MPE-seq primers (max 26bp),
# but depending on the actual primer used  in the experiment (e.g. ACT1 primer has 25bp annealing region),
# this value should be changed.
def EmptyFilter(full_path_input, full_path_output, infile, n = 33, write_every = 1000000):

	out1 = []
	out2 = []

	with gzip.open('%s/%s_r1.fastq.gz' % (full_path_input, infile)) as i1, gzip.open('%s/%s_r2.fastq.gz' % (full_path_input, infile)) as i2, gzip.open('%s/%s_r1.fastq.gz' % (full_path_output, infile) ,'wb') as o1, gzip.open('%s/%s_r2.fastq.gz' % (full_path_output, infile) , 'wb') as o2:
		i = 0
		for line1, line2 in zip_longest(i1, i2):
			
			info1 = line1.rstrip().split(b' ')[0]
			info2 = line2.rstrip().split(b' ')[0]
			

			seq1 = next(i1).rstrip()
			extra1 = next(i1).rstrip()
			quality1 = next(i1).rstrip()
			length = len(seq1)
			

			seq2 = next(i2).rstrip()
			extra2 = next(i2).rstrip()
			quality2 = next(i2).rstrip()
 
			if info1 != info2:
				print('Error: %s does not match %s' % (info1, info2))

			if length > n:
				out1.append(b'%s\n%s\n%s\n%s\n' % (info1, seq1, extra1, quality1))
				out2.append(b'%s\n%s\n%s\n%s\n' % (info2, seq2, extra2, quality2))
			
			i+=1

			if i == write_every:
				o1.write(b''.join(out1))
				o2.write(b''.join(out2))
				i = 0
				out1 = []
				out2 = []
		o1.write(b''.join(out1))
		o2.write(b''.join(out2))


main()