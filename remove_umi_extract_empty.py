import os, gzip, subprocess, sys, HTSeq
from itertools import zip_longest
from collections import defaultdict
import subprocess

def main():

	directory = sys.argv[1]
	samples = sys.argv[2:len(sys.argv)]
	for sample in samples:
		print("Removing UMI for %s"%(sample))
		ExtractMolecularCounter('%s/trimmed/'%(directory),'%s/empty/'%(directory),sample,7)

# The out put files <full_path_output/infile>_r*.fastq.gz (* is 1 or 2) differ from the input files 
# <full_path_input/infile>_r*.fastq.gzsuch that for both files the info lines only contain what is before ' ',
# plus the first <n> base umi from the r1 sequence, separated by ';'. The info between r1 and r2 are checked
# to make sure they are the same. The r1 output sequence has the first <n> bases removed.
# The file writes every <write_every> sequences processed, to adjust to the available memory.
# The default <write_every is 1000000>
# Outputs only read1 with length <= 33 + n, where n is the length of the umi. This extracts all empties.
def ExtractMolecularCounter(full_path_input, full_path_output, infile, n, write_every = 1000000):

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

			if length <= 33+n:
				out1.append(b'%s;%s\n%s\n%s\n%s\n' % (info1, seq1[:n], seq1[n:], extra1, quality1[n:]))
				out2.append(b'%s;%s\n%s\n%s\n%s\n' % (info2, seq1[:n], seq2, extra2, quality2))
			
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