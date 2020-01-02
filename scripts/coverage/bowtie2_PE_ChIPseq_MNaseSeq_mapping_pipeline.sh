#!/bin/bash

# Authors: Xiaohui Zhao (xz289@cam.ac.uk), Andrew Tock (ajt200@cam.ac.uk)
# Description: Align paired-end ChIP-seq or MNase-seq reads
#  to the Arabidopsis TAIR10 reference genome assembly and
#  filter alignments to retain uniquely aligned reads and
#  the best multiply aligned reads
# Software: bowtie2 version 2.2.9
#  samtools version 1.3
#  R version 3.3.2
#  Python version 2.7.13

# Example usage via condor submission system on hydrogen node7
# csmit -m 50G -c 24 "bash ./bowtie2_PE_ChIPseq_MNaseSeq_mapping_pipeline.sh wt_REC8_HA_Rep2_ChIP 24"

i=$1
threads=$2

gunzip --keep --force ${i}_R1.fastq.gz
gunzip --keep --force ${i}_R2.fastq.gz
# deduplicate fastq files
deduplicate.py ${i}_R1.fastq ${i}_R2.fastq
# align deduplicated reads to reference genome
bowtie2 --very-sensitive --no-discordant --no-mixed -p ${threads} -k 10 -x /projects/xiaohui_zhao/Kyuha_Spo11oligo/WT_met1_arp6_suvh_18_05_16/Fastq.files/TAIR10_chr_all -1 ${i}_R1.RmDup.fastq -2 ${i}_R2.RmDup.fastq -S ${i}_RmDup_k10_bt2.sam >> ${i}_RmDup_k10_bt2.stats 2>&1
# convert sam to bam
samtools view -@ ${threads} -bS -o ${i}_RmDup_k10_bt2.bam ${i}_RmDup_k10_bt2.sam
rm ${i}_RmDup_k10_bt2.sam
samtools view -@ ${threads} -b -hf 0x02 ${i}_RmDup_k10_bt2.bam > ${i}_RmDup_k10_bt2_mapped.bam
samtools view -@ ${threads} -h -o ${i}_RmDup_k10_bt2_mapped.sam ${i}_RmDup_k10_bt2_mapped.bam
# allow a maximum of 2 mismatches in alignment ([^0-9] matches characters not in the range of 0 to 9)
samtools view -@ ${threads} -Sh ${i}_RmDup_k10_bt2_mapped.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${i}_RmDup_k10_bt2_mapped_lowmiss.sam
rm ${i}_RmDup_k10_bt2_mapped.sam
samtools view -@ ${threads} -H ${i}_RmDup_k10_bt2_mapped_lowmiss.sam > ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_header.sam
# remove_orphans.py (a.k.a. foo.py) is by Devon Ryan,
# as posted at https://www.biostars.org/p/95929/
# and solves the problem of orphan reads in a pair that may remain after retaining
# only those that align uniquely (i.e., after filtering using grep -v "XS:i:")
samtools view -@ ${threads} -S -f 0x02 ${i}_RmDup_k10_bt2_mapped_lowmiss.sam | grep -v "XS:i:" | remove_orphans.py > ${i}_RmDup_k10_bt2_mapped_lowmiss_unique.txt
cat ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_header.sam ${i}_RmDup_k10_bt2_mapped_lowmiss_unique.txt > ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_ori.sam
rm ${i}_RmDup_k10_bt2_mapped_lowmiss_unique.txt
# bowtie2 MAPQ scores >=42 correspond to uniquely mapping reads
samtools view -@ ${threads} -h -q 42 ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_ori.sam > ${i}_RmDup_k10_bt2_mapped_lowmiss_unique.sam
rm ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_ori.sam
samtools view -@ ${threads} -bS -o ${i}_RmDup_k10_bt2_mapped_lowmiss_unique.bam ${i}_RmDup_k10_bt2_mapped_lowmiss_unique.sam
rm ${i}_RmDup_k10_bt2_mapped_lowmiss_unique.sam
samtools sort -@ ${threads} ${i}_RmDup_k10_bt2_mapped_lowmiss_unique.bam -o ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_sort.bam
samtools index ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_sort.bam

# identify and filter multiply mapping reads, and combine with uniquely mapping reads ("both")
samtools view -@ ${threads} -S -f 0x02 ${i}_RmDup_k10_bt2_mapped_lowmiss.sam | grep "XS:i:"| remove_orphans.py > ${i}_RmDup_k10_bt2_mapped_lowmiss_multi.txt
rm ${i}_RmDup_k10_bt2_mapped_lowmiss.sam
cat ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_header.sam ${i}_RmDup_k10_bt2_mapped_lowmiss_multi.txt > ${i}_RmDup_k10_bt2_mapped_lowmiss_multi.sam
rm ${i}_RmDup_k10_bt2_mapped_lowmiss_multi.txt
samtools view -@ ${threads} -h -q 10 ${i}_RmDup_k10_bt2_mapped_lowmiss_multi.sam > ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10.sam
rm ${i}_RmDup_k10_bt2_mapped_lowmiss_multi.sam
samtools view -@ ${threads} -S -f 0x02 ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10.sam | grep "XS:i:" | remove_orphans.py > ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10.txt
rm ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10.sam
wc -l ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10.txt >> ${i}_RmDup_mdim.stats 2>&1
/applications/R/R-3.3.2/bin/Rscript multi_unique_extract_pairend.R ${i}_RmDup_mdim.stats ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10.txt ${i}_RmDup_MU.RData ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique.txt
rm ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10.txt
cat ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_header.sam ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique.txt > ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique.sam
rm ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique.txt
samtools view -@ ${threads} -bS -o ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique.bam ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique.sam
rm ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique.sam ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_header.sam
samtools sort -@ ${threads} ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique.bam -o ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique_sort.bam
samtools index ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique_sort.bam
samtools merge ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_both.bam ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_sort.bam ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique_sort.bam
samtools sort -@ ${threads} ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_both.bam -o ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam
samtools index ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam

# calculate alignment summary statistics for fastq files and for each bam file generated by the pipeline
cat ${i}_R1.fastq | echo $((`wc -l`/4)) > ${i}_R1.fastq.stats
cat ${i}_R1.RmDup.fastq | echo $((`wc -l`/4)) > ${i}_R1.RmDup.fastq.stats
cat ${i}_R2.fastq | echo $((`wc -l`/4)) > ${i}_R2.fastq.stats
cat ${i}_R2.RmDup.fastq | echo $((`wc -l`/4)) > ${i}_R2.RmDup.fastq.stats
samtools view -@ ${threads} -F 0x4 ${i}_RmDup_k10_bt2.bam | cut -f 1 | sort | uniq | wc -l > ${i}_RmDup_k10_bt2.bam.stats
samtools view -@ ${threads} -F 0x4 ${i}_RmDup_k10_bt2_mapped.bam | cut -f 1 | sort | uniq | wc -l > ${i}_RmDup_k10_bt2_mapped.bam.stats
samtools view -@ ${threads} -F 0x4 ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_sort.bam | cut -f 1 | sort | uniq | wc -l > ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_sort.bam.stats
samtools view -@ ${threads} -F 0x4 ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique_sort.bam | cut -f 1 | sort | uniq | wc -l > ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique_sort.bam.stats
samtools view -@ ${threads} -F 0x4 ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam | cut -f 1 | sort | uniq | wc -l > ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam.stats
paste -d "\t" ${i}_R1.fastq.stats ${i}_R1.RmDup.fastq.stats ${i}_RmDup_k10_bt2_mapped.bam.stats ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_sort.bam.stats ${i}_RmDup_k10_bt2_mapped_lowmiss_multi_fq10_unique_sort.bam.stats ${i}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam.stats > ${i}_alignment_summary_stats.txt
sed -i '1i Total_sequenced_read_pairs\tDeduplicated\tAligned\tUniquely_aligning,_mismatches<=2\tMultiply_aligning,_mismatches<=2,_MAPQ>=10\tBoth' ${i}_alignment_summary_stats.txt 
[ -d alignment_stats ] || mkdir alignment_stats
mv ${i}*.stats alignment_stats/
mv ${i}_alignment_summary_stats.txt alignment_stats/
# delete uncompressed fastq files
rm ${i}_R1.fastq ${i}_R1.RmDup.fastq ${i}_R2.fastq ${i}_R2.RmDup.fastq
