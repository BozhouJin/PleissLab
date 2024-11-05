import os, gzip, subprocess, sys, HTSeq
from itertools import zip_longest
from collections import defaultdict
import subprocess

def main():

	directory = "/home/jbz/Desktop/PleissLab/Illumina/2.062023FragTemp/10.072024Multi"
	samples = ["471","491","511","531","472","492","512","532"]
	for sample in samples:
		print("Removing primer for %s"%(sample))
		ExtractMolecularCounter('%s'%(directory),'%s'%(directory),sample)

# Trims the first 26 bases of unmapped R1
# Modified from remove_umi_filter_empty.py 
def ExtractMolecularCounter(full_path_input, full_path_output, infile, write_every = 1000000):

	out1 = []
	out2 = []

	with open('%s/%s_Unmapped.out.mate1' % (full_path_input, infile)) as i1, gzip.open('%s/%s_r1.fastq.gz' % (full_path_output, infile) ,'wb') as o1:
		i = 0
		for line1 in i1:
			
			info1 = line1.rstrip().split(' ')[0].encode()
			

			seq1 = next(i1).rstrip().encode()
			extra1 = next(i1).rstrip().encode()
			quality1 = next(i1).rstrip().encode()


			
			out1.append(b'%s\n%s\n%s\n%s\n' % (info1, seq1[26:], extra1, quality1[26:]))

			i+=1

			if i == write_every:
				o1.write(b''.join(out1))

				i = 0
				out1 = []

		o1.write(b''.join(out1))



main()