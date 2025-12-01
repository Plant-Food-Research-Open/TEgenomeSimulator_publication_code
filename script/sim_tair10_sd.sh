#!/bin/bash

##########################
## scenario I 
## high idn 90-100 
## high sd 15-20
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_sd_scenario1
minidn=90
maxidn=100
minsd=15
maxsd=20
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_scenario1
TIME="1-00:00:00"
THREADS=1
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## scenario II
## low idn 70-80 
## high sd 15-20
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_sd_scenario2
minidn=70
maxidn=80
minsd=15
maxsd=20
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_scenario2
TIME="1-00:00:00"
THREADS=1
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## scenario III
## low idn 70-80 
## low sd 1-5
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_sd_scenario3
minidn=70
maxidn=80
minsd=1
maxsd=5
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_scenario3
TIME="1-00:00:00"
THREADS=1
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## scenario IV
## high idn 90-100 
## low sd 1-5
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_sd_scenario4
minidn=90
maxidn=100
minsd=1
maxsd=5
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_scenario4
TIME="1-00:00:00"
THREADS=1
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'