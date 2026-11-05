#!/bin/bash

##########
## run7 ##
##########
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_run7
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEgenomeSimulator
TIME="01:00:00"
THREADS=20
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -S $tetools -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

