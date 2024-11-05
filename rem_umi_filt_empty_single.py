import os, gzip, subprocess, sys, HTSeq
from itertools import zip_longest
from collections import defaultdict
import subprocess

def main():

	directory = sys.argv[1]
	samples = sys.argv[2:len(sys.argv)]
	for sample in samples:
		print("Removing UMI for %s"%(sample))
		ExtractMolecularCounter('%s/trimmed/'%(directory),'%s/noUmiNoEmpty/'%(directory),sample,7)

# 2022.8.30: This file was modified from remove_umi_filter_empty.py to do the same tasks on single-end reads
def ExtractMolecularCounter(full_path_input, full_path_output, infile, n, write_every = 1000000):

	out1 = []

	with gzip.open('%s/%s_r1.fastq.gz' % (full_path_input, infile)) as i1, gzip.open('%s/%s_r1.fastq.gz' % (full_path_output, infile) ,'wb') as o1:
		i = 0
		for line1 in i1:
			
			info1 = line1.rstrip().split(b' ')[0]
			seq1 = next(i1).rstrip()
			extra1 = next(i1).rstrip()
			quality1 = next(i1).rstrip()
			length = len(seq1)

			if length > 33+n:
				out1.append(b'%s;%s\n%s\n%s\n%s\n' % (info1, seq1[:n], seq1[n:], extra1, quality1[n:]))
			
			i+=1

			if i == write_every:
				o1.write(b''.join(out1))
				i = 0
				out1 = []

		o1.write(b''.join(out1))



main()