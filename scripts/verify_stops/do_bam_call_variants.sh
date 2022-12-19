#!/bin/bash
#$ -v SGE_FACILITIES
#$ -P unified
#$ -l h_rt=259200,h_vmem=24G,mem_free=24G,m_mem_free=8G
#$ -pe multicore 4

source ~/.bashrc
conda activate lofreq


alltasks=('Deer_1271_S63' 'Deer_2070_S44' 'Deer_1839_S75')

bwa mem -t 4 ~/SARS-CoV2_reference_NC_045512.fasta ${alltasks[$SGE_TASK_ID-1]}_L001_R1_001.fastq.gz ${alltasks[$SGE_TASK_ID-1]}_L001_R2_001.fastq.gz > ${alltasks[$SGE_TASK_ID-1]}.sam

samtools view -b -S -o ${alltasks[$SGE_TASK_ID-1]}.bam ${alltasks[$SGE_TASK_ID-1]}.sam
samtools sort -o ${alltasks[$SGE_TASK_ID-1]}_sorted.bam ${alltasks[$SGE_TASK_ID-1]}.bam 
samtools index ${alltasks[$SGE_TASK_ID-1]}_sorted.bam

LOFREQ=lofreq
REF=~/SARS-CoV2_reference_NC_045512.fasta
PREFIX=${alltasks[$SGE_TASK_ID-1]}
BAM=${PREFIX}_sorted.bam
VITERBAM=${PREFIX}_viterbi.bam
SORTED=${PREFIX}_viterbi_sorted.bam
IQBAM=${PREFIX}_viterbi_iq.bam
VCF=${PREFIX}.vcf
SNVS=${PREFIX}.dp4.SNVs.vcf
INDELS=${PREFIX}.dp10.indels.vcf

$LOFREQ viterbi -f $REF -o $VITERBAM $BAM
samtools sort -o $SORTED $VITERBAM
$LOFREQ indelqual -f $REF -o $IQBAM --dindel ${SORTED}
samtools index $IQBAM

#call variants using  min depth 2
$LOFREQ call -f $REF -o $VCF  -C 2 --call-indels --no-default-filter --force-overwrite --use-orphan $IQBAM

#call indels
$LOFREQ filter  -Q 20 -K 20 --no-defaults  -v 4 -V 0 -a 0.500001 -A 0 --only-snvs -i $VCF -o $SNVS

#do indels with DP >= 20 and AF >= 0.5
$LOFREQ filter  -Q 20 -K 20 --no-defaults  -v 10 -V 0 -a 0.50000 -A 0 --only-indels -i  $VCF -o $INDELS
