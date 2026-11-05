#!/bin/bash
# GITHUB_REPO=cloned repository of "TEgenomeSimulator"
# Navigate to the cloned repository of "TEgenomeSimulator_publication_code"
ml apptainer
ml samtools/1.20
#TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # remove this line before git push
cd $WRK
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEGS_${prefix}
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"  # remove this line before git push
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
# Submitted batch job 10737364

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

