#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --job-name=hic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=200GB

module load bowtie/2.2.9
module load R/3.5.0
###############################
# Build bowtie2 indices
###############################
bowtie2-build karyon_assignment_98_210/210V2AP2_karyon_ALL.fa bowtie_indices/210V2AP2_karyon_ALL
###############################
# Create digested reference genome
###############################
# The argument ‘–re1’ specifies the restriction enzyme used to digest the genome 
# (a caret symbol ‘^’ is used to denote the restriction enzyme cut site, and a 
# comma separates the DNA sequence from the restriction enzyme name). 
###############################
./hicup_v0.7.1/hicup_digester --genome 210V2AP2_karyon_ALL --re1 ^GATC,Sau3A *.fa --outdir digester_out
###############################
module load samtools
###############################
mkdir hic_analysis
./hicup_v0.7.1/hicup --config hicup.conf --threads 8 --outdir hic_analysis
###############################
