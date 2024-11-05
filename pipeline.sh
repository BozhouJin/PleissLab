#!/bin/bash

# Set the home directory and the sample names and the raw names
# The home directory is upper to the raw data directory
# The sample names will be used to generate the resulting data files
# The raw names are the names used to distinguish the raw sequence data files
samples=( "41" "42" "AILS1" "AILS2")
raw_names=( "41" "42" "AILS1" "AILS2")
home_dir="/home/jbz/Desktop/PleissLab/Illumina/5.112024AILSoffTarget"


# If not already, create the necessary directories. Use by removeing the comment block (the "<< pair"), same for below
# 20240424: added a folder for empties
<<mkdir
mkdir "${home_dir}/fastqc" "${home_dir}/aligned" "${home_dir}/noUmiNoEmpty" "${home_dir}/trimmed" "${home_dir}/empty"
mkdir



# Quality check. 
<<Qcheck
echo "Generating fastqc files for each datafile"
for file in `ls "${home_dir}/raw"`
do
	fastqc "${home_dir}/raw/${file}" -o "${home_dir}/fastqc"
done
Qcheck



# Adapter trimming. This removes the adapter on the 3' side of each read
# due to short insert and long reads.
<<TrimA
echo "Trimming adapters"

for ((i=0;i<${#samples[@]};i++))
do
	fastp -i "${home_dir}/raw/"*_${raw_names[$i]}_*R1* -I "${home_dir}/raw/"*_${raw_names[$i]}_*R2* -o "${home_dir}/trimmed/${samples[$i]}_r1.fastq.gz" -O "${home_dir}/trimmed/${samples[$i]}_r2.fastq.gz" --disable_trim_poly_g --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT  --unpaired1 "${home_dir}/trimmed/${samples[$i]}_r1_unpaired.fastq.gz" --unpaired2 "${home_dir}/trimmed/${samples[$i]}_r2_unpaired.fastq.gz" --html "${home_dir}/trimmed/${samples[$i]}.html" --json "${home_dir}/trimmed/${samples[$i]}.json" --thread 1
done
TrimA
# 2022.8.31 edit: Single-end read version
<<TrimASingle
echo "Trimming adapters"

for ((i=0;i<${#samples[@]};i++))
do
	fastp -i "${home_dir}/raw/"*_${raw_names[$i]}_*R1* -o "${home_dir}/trimmed/${samples[$i]}_r1.fastq.gz" --disable_trim_poly_g --adapter_sequence CTGTCTCTTATACACATCT --html "${home_dir}/trimmed/${samples[$i]}.html" --json "${home_dir}/trimmed/${samples[$i]}.json" --thread 1
done
TrimASingle






# Remove UMI. This is 7bp of UMI in the MPE-seq primer between the annealing sequence and
# the illumina adapter. Originally used to distinguish between individual RNA molecules,
# but was not used
# 2022.7.5 edit: this function now also filters empties (read pairs with r1 <= 33bp)
# Which probably also removes some short reads
<<Remumiempt
python3 remove_umi_filter_empty.py ${home_dir} ${samples[@]}
Remumiempt
# 2022.8.31 edit: single-end read version
<<RUFESingle
python3 rem_umi_filt_empty_single.py ${home_dir} ${samples[@]}
RUFESingle
# Empty filtering without removing UMI. Basically the same as above, but since SPE-seq primers
# do not have UMI, this is specifically designed for them.
<<Empfilt
python3 filter_empty.py ${home_dir} ${samples[@]}
Empfilt
# 20240424: Extracting empties but not extended reads. Comes with a remove umi version and one that does not
<<RemoveUmiExtractEmpty
python3 remove_umi_extract_empty.py ${home_dir} ${samples[@]}
RemoveUmiExtractEmpty
<<ExtractEmpty
python3 extract_empty.py ${home_dir} ${samples[@]}
ExtractEmpty
# The single read version:
<<RMumiExtractEmptySingle
python3 remove_umi_extract_empty_single.py ${home_dir} ${samples[@]}
RMumiExtractEmptySingle





# Generate the genome index. This only needs to be done once in theory, but needs to be updated
# when, e.g., read length changes. The original readlength 151-7umi-1 = 143, so change
# sjdbOverhang to 143 according to the manual
<<Index
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /home/jbz/Desktop/PleissLab/pombe_cerevisiae_index/index --genomeFastaFiles /home/jbz/Desktop/PleissLab/pombe_cerevisiae_index/concat_genome.fa --sjdbGTFfile /home/jbz/Desktop/PleissLab/pombe_cerevisiae_index/concat_features.gff --sjdbGTFtagExonParentTranscript Parent  --sjdbOverhang 150
Index




# Sequence alignment with STAR
<<Align
for sample in ${samples[@]}
do
	echo "Aligning ${sample}"
	STAR --genomeDir /home/jbz/Desktop/PleissLab/pombe_cerevisiae_index/index --readFilesIn ${home_dir}/noUmiNoEmpty/${sample}_r1.fastq.gz ${home_dir}/noUmiNoEmpty/${sample}_r2.fastq.gz  --readFilesCommand zcat --peOverlapNbasesMin 5 --peOverlapMMp 0.1 --outFilterMultimapNmax 1 --alignIntronMin 10 --alignIntronMax 1100 --outSAMattributes All --runThreadN 4 --outSAMunmapped Within KeepPairs --alignSJoverhangMin 3 --alignSplicedMateMapLmin 3 --alignMatesGapMax 3000 --outFilterMismatchNmax 999 --alignEndsType EndToEnd --outFileNamePrefix ${home_dir}/aligned/${sample}_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
	samtools index ${home_dir}/aligned/${sample}_Aligned.sortedByCoord.out.bam
done
Align
# 2022.8.30 edit: single-end version
<<AlignSingle
for sample in ${samples[@]}
do
	echo "Aligning ${sample}"
	STAR --genomeDir /home/jbz/Desktop/PleissLab/pombe_cerevisiae_index/index --readFilesIn ${home_dir}/noUmiNoEmpty/${sample}_r1.fastq.gz --readFilesCommand zcat --peOverlapNbasesMin 5 --peOverlapMMp 0.1 --outFilterMultimapNmax 1 --alignIntronMin 10 --alignIntronMax 1100 --outSAMattributes All --runThreadN 4 --outSAMunmapped Within KeepPairs --alignSJoverhangMin 3 --alignSplicedMateMapLmin 3 --alignMatesGapMax 3000 --outFilterMismatchNmax 999 --alignEndsType EndToEnd --outFileNamePrefix ${home_dir}/alignedSingle/${sample}_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
	samtools index ${home_dir}/alignedSingle/${sample}_Aligned.sortedByCoord.out.bam
done
AlignSingle
# 20230504 edit: this is alignment using R1 only, and allowing multimapping.
<<AlignSingle_multiAllowed
for sample in ${samples[@]}
do
	echo "Aligning ${sample}"
	STAR --genomeDir /home/jbz/Desktop/PleissLab/pombe_cerevisiae_index/index --readFilesIn ${home_dir}/noUmiNoEmpty/${sample}_r1.fastq.gz --readFilesCommand zcat --peOverlapNbasesMin 5 --peOverlapMMp 0.1 --outFilterMultimapNmax 9999 --alignIntronMin 10 --alignIntronMax 1100 --outSAMattributes All --runThreadN 4 --outSAMunmapped Within KeepPairs --alignSJoverhangMin 3 --alignSplicedMateMapLmin 3 --alignMatesGapMax 3000 --outFilterMismatchNmax 999 --alignEndsType EndToEnd --outFileNamePrefix ${home_dir}/alignedSingle_multiAllowed/${sample}_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
	samtools index ${home_dir}/alignedSingle_multiAllowed/${sample}_Aligned.sortedByCoord.out.bam
done
AlignSingle_multiAllowed



# Feature counting
<<FeatC
python3 feature_count.py ${home_dir} ${samples[@]}
FeatC
# 2022.8.30 edit: single-end version
<<FeatSingle
python3 feature_single.py ${home_dir} ${samples[@]}
FeatSingle

# Empty counting
#<<Empty
python3 empty_count.py ${home_dir} ${samples[@]}
#Empty

# Empty features counting
<<EmptyFeatures
python3 empty_feature_count.py ${home_dir} ${samples[@]}
EmptyFeatures