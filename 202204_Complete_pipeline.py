import os, gzip, subprocess, HTSeq, sys
from itertools import izip_longest
from collections import defaultdict
import subprocess

def main():

	samples = []

	home_directory = '/media/8T/srd96/Karla/202204_Timepoint6_Screen'

	#makes list of samples in home directory
	for file in os.listdir('%s/Raw' % (home_directory)):
		if '_R1' in file:
			file = file.split('_R')[0]
			samples.append(file)
	print samples


	for sample in samples:
		print sample
		qualityCheck('%s/Raw' % (home_directory), sample)
		print 'Fastqc Generated'
		trimAdapters('%s/Raw'% (home_directory), '%s/Trim'% (home_directory), sample)
		print 'Trimmed'
		ExtractMolecularCounter('%s/Trim' % (home_directory), '%s/NoUMI' % (home_directory), sample, 7)
		print 'UMI Extracted'
		Alignment(home_directory, sample)
		print 'Aligned'

	#start Feature counting
	# Read in intron ranges
	targets, introns, orientation = readTargetFeatures('/home/srd96/Genomes/concat_intron_ranges.bed')
	bp = readBranchPoints('/home/srd96/Genomes/Branch_window.txt')

	
	intermediate_counts = defaultdict(int)

	print "Feature Counting..."

	for sample in samples:

		feature_count('%s/Align' % (home_directory), sample, targets, introns, orientation, bp, intermediate_counts)

	with open('%s/Count/Total_Counts_20220426.txt' % (home_directory), 'w') as outfile:
		outfile.write('Sample\tIntron\tMature\tPre_First_Step\tLariat\tUnknown\tPolyG\tBad_Read2\tPremature_Total\n')
		for target in introns:
			for sample in samples:
				outfile.write('%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % (sample, target, intermediate_counts[(sample, target, 'Mature')], intermediate_counts[(sample, target, 'Pre_First_Step')], intermediate_counts[(sample, target, 'Lariat')], intermediate_counts[(sample, target, 'Unknown')], intermediate_counts[(sample, target, 'PolyG')], intermediate_counts[(sample, target, 'Bad_Read2')], intermediate_counts[(sample, target, 'Premature_Total')]))
	
	#start empty counting

	print "Empty Counting..."

	counts = defaultdict(float)

	for sample in samples:
		EmptyCounting(sample, counts, home_directory)

	with open('%s/Count/EmptyCounts_20220426.txt' % (home_directory), 'w') as outfile:
		outfile.write('Sample\tTotal\tEmpty\tPercentEmpty\n')
		for sample in samples:
			outfile.write('%s\t%f\t%f\t%f\n' % (sample, counts[(sample, 'total')], counts[(sample, 'empty')], (counts[(sample, 'empty')]/counts[(sample, 'total')])*100))


def qualityCheck(full_path, infile):
	cmd = 'fastqc %s/%s_R1.fastq.gz' % (full_path, infile)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()

	cmd = 'fastqc %s/%s_R2.fastq.gz' % (full_path, infile)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()

def trimAdapters(full_path_raw, full_path_output, infile):
	cmd = 'fastp -i %s/%s_R1.fastq.gz -I %s/%s_R2.fastq.gz -o %s/%s_Trimmed_R1_P.fastq.gz -O %s/%s_Trimmed_R2_P.fastq.gz --disable_trim_poly_g --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT  --unpaired1 %s/%s_unpaired_Trim.fastq.gz --unpaired2 %s/%s_unpaired_Trim.fastq.gz --html %s/%s.html --json %s/%s.json --thread 2' % (full_path_raw, infile, full_path_raw, infile, full_path_output, infile, full_path_output, infile, full_path_output, infile, full_path_output, infile, full_path_output, infile, full_path_output, infile)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()

def ExtractMolecularCounter(full_path_input, full_path_output, infile, n):

	out1 = []
	out2 = []

	with gzip.open('%s/%s_Trimmed_R1_P.fastq.gz' % (full_path_input, infile)) as i1, gzip.open('%s/%s_Trimmed_R2_P.fastq.gz' % (full_path_input, infile)) as i2:
		for line1, line2 in izip_longest(i1, i2):
			
			info1 = line1.rstrip().split(' ')[0]
			info2 = line2.rstrip().split(' ')[0]
			

			seq1 = next(i1).rstrip()
			extra1 = next(i1).rstrip()
			quality1 = next(i1).rstrip()
			

			seq2 = next(i2).rstrip()
			extra2 = next(i2).rstrip()
			quality2 = next(i2).rstrip()
 
			if info1 != info2:
				print 'Error: %s does not match %s' % (info1, info2)


			out1.append('%s;%s\n%s\n%s\n%s\n' % (info1, seq1[:n], seq1[n:], extra1, quality1[n:]))
			out2.append('%s;%s\n%s\n%s\n%s\n' % (info2, seq1[:n], seq2, extra2, quality2))
				
	with gzip.open('%s/%s_R1_NoUMI.fastq.gz' % (full_path_output, infile) ,'wb') as out:
		out.write(''.join(out1))

	with gzip.open('%s/%s_R2_NoUMI.fastq.gz' % (full_path_output, infile) , 'wb') as out:
		out.write(''.join(out2))


def Alignment(path, infile):
	
	cmd = 'STAR --genomeDir /home/srd96/Genomes/STAR_genome --readFilesIn %s/NoUMI/%s_R1_NoUMI.fastq.gz %s/NoUMI/%s_R2_NoUMI.fastq.gz  --readFilesCommand zcat --peOverlapNbasesMin 5 --peOverlapMMp 0.1 --outFilterMultimapNmax 1 --alignIntronMin 10 --alignIntronMax 1100 --outSAMattributes All --runThreadN 4 --outSAMunmapped Within KeepPairs --alignSJoverhangMin 3 --alignSplicedMateMapLmin 3 --alignMatesGapMax 3000 --outFilterMismatchNmax 999 --alignEndsType EndToEnd --outFileNamePrefix %s/Align/%s_ --outSAMtype BAM SortedByCoordinate' % (path, infile, path, infile, path, infile)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()

	cmd = 'samtools index %s/Align/%s_Aligned.sortedByCoord.out.bam' % (path, infile)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()

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

	#How does this pairing know which read is the original read1? Because of the initial sorting in alignment?
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
		
				#If read1 maps with >= 1 junctions, count it as mature?
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
					elif read2 != None and read2.aligned and read2.read == 'GGGGGGGGGGGGGGGGG':
						intermediate_counts[(sample, target, 'PolyG')] += 1
					else:
						intermediate_counts[(sample, target, 'Bad_Read2')] += 1
	
def EmptyCounting(sample, counts, home_directory):

	with gzip.open("%s/NoUMI/%s_R1_NoUMI.fastq.gz" % (home_directory, sample)) as infile:
		
		

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


main()