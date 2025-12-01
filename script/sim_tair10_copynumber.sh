#!/bin/bash

######################
## copy number 1-10 ##
######################
# Navigate to the cloned repository of TEgenomeSimulator
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_cn_1_10
min=1
max=10
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#######################
## copy number 5-100 ##
#######################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_cn_5_100
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#######################
## copy number 5-500 ##
#######################
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_cn_5_500
min=5
max=500
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

########################
## copy number 5-1000 ##
########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_cn_5_1000
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

########################
## copy number 5-2000 ##
########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_cn_5_2000
min=5
max=2000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'