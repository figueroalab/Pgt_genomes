#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --job-name=salsa
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=200GB

export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
#---------------------------------------
#---------------------------------------
module load bowtie/1.2.0
module load samtools/1.9.0
module load bedtools
module load bwa
module load picard
module load salsa/2.2
#####################################################
# SALSA was run three times using diferent combinations
# of the A-contigs, B-contigs and ?-contigs
#####################################################
#fasta='A_Q_contigs'
#fastafile='A_Q_contigs.fasta'
#faifile='A_Q_contigs.fasta.fai'

#fasta='B_Q_contigs'
#fastafile='B_Q_contigs.fasta'
#faifile='B_Q_contigs.fasta.fai'

fasta='A_B_Q_contigs'
fastafile='A_B_Q_contigs.fasta'
faifile='A_B_Q_contigs.fasta.fai'

path='/datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/'
base='/flush1/username/Pgt_HiC_'
perlfilter='/home/username/Pgt_210_HiC_Salsa_Scaffolding/filter_five_end.pl'
perlcombiner='/home/username/Pgt_210_HiC_Salsa_Scaffolding/two_read_bam_combiner.pl'
#####################################################
#####################################################
cat /datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/A_contigs.fasta /datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/B_contigs.fasta /datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/Question_mark_contigs.fasta > /datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/A_B_Q_contigs.fasta
cat /datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/A_contigs.fasta /datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/Question_mark_contigs.fasta > /datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/A_Q_contigs.fasta
cat /datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/B_contigs.fasta /datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/Question_mark_contigs.fasta > /datastore/username/Pgt_21_0_Data/karyon_assignment_98_210/B_Q_contigs.fasta
#####################################################
echo $base$fasta
mkdir $base$fasta
cd $base$fasta
#####################################################
bwa index -a bwtsw -p $fasta $path$fastafile

bwa mem -t4 $fasta /datastore/username/Pgt_21_0_Data/HiC_Data/PgtHIC_S1_L001_R1_001.fastq.gz | samtools view -Sb - > PgtHIC_S1_L001_R1_001.bam
bwa mem -t4 $fasta /datastore/username/Pgt_21_0_Data/HiC_Data/PgtHIC_S1_L001_R2_001.fastq.gz | samtools view -Sb - > PgtHIC_S1_L001_R2_001.bam
bwa mem -t4 $fasta /datastore/username/Pgt_21_0_Data/HiC_Data/PgtHIC_S1_L002_R1_001.fastq.gz | samtools view -Sb - > PgtHIC_S1_L002_R1_001.bam
bwa mem -t4 $fasta /datastore/username/Pgt_21_0_Data/HiC_Data/PgtHIC_S1_L002_R2_001.fastq.gz | samtools view -Sb - > PgtHIC_S1_L002_R2_001.bam
bwa mem -t4 $fasta /datastore/username/Pgt_21_0_Data/HiC_Data/PgtHIC_S1_L003_R1_001.fastq.gz | samtools view -Sb - > PgtHIC_S1_L003_R1_001.bam
bwa mem -t4 $fasta /datastore/username/Pgt_21_0_Data/HiC_Data/PgtHIC_S1_L003_R2_001.fastq.gz | samtools view -Sb - > PgtHIC_S1_L003_R2_001.bam
bwa mem -t4 $fasta /datastore/username/Pgt_21_0_Data/HiC_Data/PgtHIC_S1_L004_R1_001.fastq.gz | samtools view -Sb - > PgtHIC_S1_L004_R1_001.bam
bwa mem -t4 $fasta /datastore/username/Pgt_21_0_Data/HiC_Data/PgtHIC_S1_L004_R2_001.fastq.gz | samtools view -Sb - > PgtHIC_S1_L004_R2_001.bam
#####################################################
samtools faidx $path$fastafile

mkdir filtered
samtools view -h PgtHIC_S1_L001_R1_001.bam | perl $perlfilter | samtools view -Sb - > filtered/PgtHIC_S1_L001_R1_001.bam
samtools view -h PgtHIC_S1_L001_R2_001.bam | perl $perlfilter | samtools view -Sb - > filtered/PgtHIC_S1_L001_R2_001.bam
samtools view -h PgtHIC_S1_L002_R1_001.bam | perl $perlfilter | samtools view -Sb - > filtered/PgtHIC_S1_L002_R1_001.bam
samtools view -h PgtHIC_S1_L002_R2_001.bam | perl $perlfilter | samtools view -Sb - > filtered/PgtHIC_S1_L002_R2_001.bam
samtools view -h PgtHIC_S1_L003_R1_001.bam | perl $perlfilter | samtools view -Sb - > filtered/PgtHIC_S1_L003_R1_001.bam
samtools view -h PgtHIC_S1_L003_R2_001.bam | perl $perlfilter | samtools view -Sb - > filtered/PgtHIC_S1_L003_R2_001.bam
samtools view -h PgtHIC_S1_L004_R1_001.bam | perl $perlfilter | samtools view -Sb - > filtered/PgtHIC_S1_L004_R1_001.bam
samtools view -h PgtHIC_S1_L004_R2_001.bam | perl $perlfilter | samtools view -Sb - > filtered/PgtHIC_S1_L004_R2_001.bam

perl $perlcombiner filtered/PgtHIC_S1_L001_R1_001.bam filtered/PgtHIC_S1_L001_R2_001.bam samtools 10 | samtools view -bS -t $path$faifile - | samtools sort -o L001.bam -
perl $perlcombiner filtered/PgtHIC_S1_L002_R1_001.bam filtered/PgtHIC_S1_L002_R2_001.bam samtools 10 | samtools view -bS -t $path$faifile - | samtools sort -o L002.bam -
perl $perlcombiner filtered/PgtHIC_S1_L003_R1_001.bam filtered/PgtHIC_S1_L003_R2_001.bam samtools 10 | samtools view -bS -t $path$faifile - | samtools sort -o L003.bam -
perl $perlcombiner filtered/PgtHIC_S1_L004_R1_001.bam filtered/PgtHIC_S1_L004_R2_001.bam samtools 10 | samtools view -bS -t $path$faifile - | samtools sort -o L004.bam -
########
mkdir paired
picard AddOrReplaceReadGroups INPUT=L001.bam OUTPUT=paired/L001.bam ID=L001 LB=L001 SM='Pgt' PL=ILLUMINA PU=none
picard AddOrReplaceReadGroups INPUT=L002.bam OUTPUT=paired/L002.bam ID=L002 LB=L002 SM='Pgt' PL=ILLUMINA PU=none
picard AddOrReplaceReadGroups INPUT=L003.bam OUTPUT=paired/L003.bam ID=L003 LB=L003 SM='Pgt' PL=ILLUMINA PU=none
picard AddOrReplaceReadGroups INPUT=L004.bam OUTPUT=paired/L004.bam ID=L004 LB=L004 SM='Pgt' PL=ILLUMINA PU=none

INPUTS_TECH_REPS=('INPUT=L001.bam' 'INPUT=L002.bam' 'INPUT=L003.bam' 'INPUT=L004.bam')
picard MergeSamFiles $INPUTS_TECH_REPS OUTPUT=temp.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

picard MarkDuplicates INPUT=temp.bam OUTPUT=alignment.bam METRICS_FILE=metrics.txt ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index alignment.bam
bamToBed -i alignment.bam > alignment.bed
sort -k 4 alignment.bed > tmp && mv tmp alignment.bed

python $SALSA_HOME/run_pipeline.py -a $path$fastafile -l $path$faifile -b alignment.bed -e GATC -o scaffolds
