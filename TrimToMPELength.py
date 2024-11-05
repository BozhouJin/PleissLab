import os, gzip
from itertools import izip_longest

def main():

	samples = []
	
	home_directory = '/media/8T/srd96/Karla/202203_SpikeIn'

	for file in os.listdir('%s/Raw' % (home_directory)):
		if '_R1' in file:
			file = file.split('_R')[0]
			samples.append(file)
	print samples

	for sample in samples:

		out1 = []
		out2 = []
	
	
		with gzip.open('%s/Raw/%s_R1.fastq.gz' % (home_directory, sample)) as i1, gzip.open('%s/Raw/%s_R2.fastq.gz' % (home_directory, sample)) as i2:
			for line1, line2 in izip_longest(i1, i2):
				
				info1 = line1.rstrip().split(' ')[0]
				info2 = line2.rstrip().split(' ')[0]
				
	
				seq1 = next(i1).rstrip()[:62]
				extra1 = next(i1).rstrip()
				quality1 = next(i1).rstrip()[:62]
				
	
				seq2 = next(i2).rstrip()[:17]
				extra2 = next(i2).rstrip()
				quality2 = next(i2).rstrip()[:17]
 	
				if info1 != info2:
					print 'Error: %s does not match %s' % (info1, info2)
	
	
				out1.append('%s\n%s\n%s\n%s\n' % (info1, seq1, extra1, quality1))
				out2.append('%s\n%s\n%s\n%s\n' % (info2, seq2, extra2, quality2))
					
		with gzip.open('%s/%s_R1_MPE_Length.fastq.gz' % ('/media/8T/srd96/Karla/202203_SpikeIn/Raw/MPE_Length', sample) ,'wb') as out:
			out.write(''.join(out1))
	
		with gzip.open('%s/%s_R2_MPE_Length.fastq.gz' % ('/media/8T/srd96/Karla/202203_SpikeIn/Raw/MPE_Length', sample) , 'wb') as out:
			out.write(''.join(out2))

main()