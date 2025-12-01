#!/bin/bash

##########################
### alpha 0.5 beta 0.5 ###
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_itg_1
min=5
max=100
alpha=0.5
beta=0.5
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_itg_1
TIME="00:10:00"
THREADS=1
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -a $alpha -b $beta -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

###########################
### alpha 0.5 beta 0.75 ###
###########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_itg_2
min=5
max=100
alpha=0.5
beta=0.75
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_itg_2
TIME="00:10:00"
THREADS=1
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -a $alpha -b $beta -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

########################
### alpha 0.5 beta 1 ###
########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_itg_3
min=5
max=100
alpha=0.5
beta=1
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_itg_3
TIME="00:10:00"
THREADS=1
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -a $alpha -b $beta -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

###########################
### alpha 0.75 beta 0.5 ###
###########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_itg_4
min=5
max=100
alpha=0.75
beta=0.5
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_itg_4
TIME="00:10:00"
THREADS=1
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -a $alpha -b $beta -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


########################
### alpha 1 beta 0.5 ###
########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_itg_5
min=5
max=100
alpha=1
beta=0.5
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_itg_5
TIME="00:10:00"
THREADS=1
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -a $alpha -b $beta -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'