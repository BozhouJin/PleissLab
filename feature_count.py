import os, gzip, subprocess, sys, HTSeq
from itertools import zip_longest
from collections import defaultdict
import subprocess


def readTargetFeatures(interval):
	intron_set = set()
	orientation = {}
	targets = HTSeq.GenomicArrayOfSets('auto', stranded=False)
	for line in open(interval):
		fields = line.rstrip().split('\t')
		iv = HTSeq.GenomicInterval(fields[0], int(fields[1])-1, int(fields[2]), fields[5])
		targets[iv] += fields[3]
		intron_set.add(fields[3])
		orientation[fields[3]] = fields[5]
	return targets, intron_set, orientation

def readBranchPoints(bp_infile):
	bp = {}
	for line in open(bp_infile):
		cur = line.rstrip().split('\t')
		gene = cur[0].split(';')[0]
		rank = cur[0].split(';')[1]
		bp['%s;%s' % (gene, rank)] = (int(cur[2]), int(cur[3]))
	return bp

def feature_count(fullpath, sample, targets, intron_set, orientation, bp, intermediate_counts):
	infile = '%s/%s_Aligned.sortedByCoord.out.bam' % (fullpath, sample)

	SI_counts = defaultdict(int)
	junction_counts = defaultdict(int)
	
	#intermediate_counts = defaultdict(int)

	for read1, read2 in HTSeq.pair_SAM_alignments_with_buffer(HTSeq.BAM_Reader(infile), max_buffer_size = 6000000):
	

		#if read1 != None and read1.aligned and read1.aQual > 5 and read2 != None and read2.aligned and read2.aQual > 5:
		if read1 != None and read1.aligned and read1.aQual > 5:
			'''if len(read1.read) > 33:'''
			introns = set()
			junctions = set()


			for i, cigop in enumerate(read1.cigar):
				
				if cigop.type == 'M':
					for iv, val in targets[cigop.ref_iv].steps():
						introns |= val
				elif cigop.type == 'N':
					if read1.cigar[i-1].type =='M' and read1.cigar[i-1].size > 3 and read1.cigar[i+1].type =='M' and read1.cigar[i+1].size > 3:
						for iv, val in targets[cigop.ref_iv].steps():
							junctions |= val


			candidate_genes = set()
			for i in introns:
				candidate_genes.add(i.split(';')[0])
			for i in junctions:
				candidate_genes.add(i.split(';')[0])
		

			if len(candidate_genes) == 1:
		
				if len(junctions) >= 1:
					ranks = []
					for i in junctions:
						ranks.append(int(i.split(';')[1]))
					target = '%s;%d' % (list(junctions)[0].split(';')[0], max(ranks))
					intermediate_counts[(sample, target, 'Mature')] += 1

				elif len(junctions) == 0:
					ranks = []
					for i in introns:
						ranks.append(int(i.split(';')[1]))
					target = '%s;%d' % (list(candidate_genes)[0], max(ranks))
					intermediate_counts[(sample, target, 'Premature_Total')] += 1

					if read2 != None and read2.aligned and read2.aQual > 5:
	
						if target in bp:
								
							ori = orientation[target]
	
							if ori == '+':
								if read2.iv.start > bp[target][1]:
									intermediate_counts[(sample, target, 'Unknown')] += 1
								elif read2.iv.start > bp[target][0]:
									intermediate_counts[(sample, target, 'Lariat')] += 1
								else:
									intermediate_counts[(sample, target, 'Pre_First_Step')] += 1
							else:
								if read2.iv.end < bp[target][0]:
									intermediate_counts[(sample, target, 'Unknown')] += 1
								elif read2.iv.end < bp[target][1]:
									intermediate_counts[(sample, target, 'Lariat')] += 1
								else:
									intermediate_counts[(sample, target, 'Pre_First_Step')] += 1
						else:
							intermediate_counts[(sample,target,'Branch_Point_not_Inlcuded')] += 1
					elif read2 != None and read2.aligned and read2.read == 'GGGGGGGGGGGGGGGGG':
						intermediate_counts[(sample, target, 'PolyG')] += 1
					else:
						intermediate_counts[(sample, target, 'Bad_Read2')] += 1




def main():

	#start Feature counting
	# Read in intron ranges
	targets, introns, orientation = readTargetFeatures('/home/jbz/Desktop/PleissLab/pombe_cerevisiae_index/concat_intron_ranges.bed')
	bp = readBranchPoints('/home/jbz/Desktop/PleissLab/pombe_cerevisiae_index/Branch_window.txt')

	intermediate_counts = defaultdict(int)

	print("Feature Counting...")

	directory = sys.argv[1]
	samples = sys.argv[2:len(sys.argv)]

	for sample in samples:
		print("Counting %s"%(sample))
		feature_count('%s/aligned'%(directory), sample, targets, introns, orientation, bp, intermediate_counts)

	with open('%s/Total_Counts.txt'%(directory), 'wt') as outfile:
		outfile.write('Sample\tIntron\tMature\tPre_First_Step\tLariat\tUnknown\tPolyG\tBad_Read2\tBranch_Point_not_Included\tPremature_Total\n')
		for target in introns:
			for sample in samples:
				outfile.write('%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % (sample, target, intermediate_counts[(sample, target, 'Mature')], intermediate_counts[(sample, target, 'Pre_First_Step')], intermediate_counts[(sample, target, 'Lariat')], intermediate_counts[(sample, target, 'Unknown')], intermediate_counts[(sample, target, 'PolyG')], intermediate_counts[(sample, target, 'Bad_Read2')],intermediate_counts[(sample, target, 'Branch_Point_not_Inlcuded')], intermediate_counts[(sample, target, 'Premature_Total')]))
	



main()