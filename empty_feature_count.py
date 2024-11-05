# This file counts the empty reads mapped to each gene
import os, gzip, subprocess, sys, HTSeq
from itertools import zip_longest
from collections import defaultdict
import subprocess


def readTargetFeatures(interval):
	orientation = {}
	intron_set = set()
	targets = HTSeq.GenomicArrayOfSets('auto',stranded=False)
	for line in open(interval):
		fields = line.rstrip().split('\t')
		iv = HTSeq.GenomicInterval(fields[0], int(fields[1])-11, int(fields[2])+10, fields[5])
		targets[iv] += fields[3]
		intron_set.add(fields[3])
		orientation[fields[3]] = fields[5]
	return targets, intron_set, orientation


def feature_count(fullpath, sample, targets, orientation,counts):
	infile = '%s/%s_Aligned.sortedByCoord.out.bam' % (fullpath, sample)

	for read1 in HTSeq.BAM_Reader(infile):
		# Did not include mapping quality check, which should allow multimapping
		# of primers that target paralogs.
		if read1 != None and read1.aligned:
			introns = set()
			for i, cigop in enumerate(read1.cigar):
				if cigop.type == 'M':
					for iv, val in targets[cigop.ref_iv].steps():
						introns |= val
			for intron in introns:
				counts[(sample, intron)] += 1



def main():

	#start Feature counting
	# Read in intron ranges
	targets, intron_set, orientation = readTargetFeatures('/home/jbz/Desktop/PleissLab/pombe_cerevisiae_index/concat_intron_ranges.bed')

	print("Feature Counting...")
	counts = defaultdict(int)


	directory = sys.argv[1]
	samples = sys.argv[2:len(sys.argv)]

	for sample in samples:
		print("Counting %s"%(sample))
		feature_count('%s/empty_mapped'%(directory), sample, targets, orientation,counts)


	with open('%s/Empty_features.txt'%(directory), 'wt') as outfile:
		outfile.write('Sample\tGene\tEmpty\n')
		for intron in intron_set:
			for sample in samples:
				if counts[(sample, intron)] > 0:
					outfile.write('%s\t%s\t%d\n' % (sample, intron, counts[(sample, intron)]))



main()