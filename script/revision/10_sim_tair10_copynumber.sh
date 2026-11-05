#!/bin/bash

######################
## copy number 1-10 ##
######################
# GITHUB_REPO=cloned repository of "TEgenomeSimulator"
# Navigate to the cloned repository of "TEgenomeSimulator_publication_code"
ml apptainer
ml samtools/1.20
#TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # remove this line before git push
cd $WRK
#TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_cn_1_10
min=1
max=10
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEGS_${prefix}
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"  # remove this line before git push
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"
# Submitted batch job 10737369

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#######################
## copy number 5-100 ##
#######################
prefix=tair10_cn_5_100
min=5
max=100

JOB=TEGS_${prefix}
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"
# Submitted batch job 10737370

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#######################
## copy number 5-500 ##
#######################
prefix=tair10_cn_5_500
min=5
max=500

JOB=TEGS_${prefix}
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"
# Submitted batch job 10737371

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

########################
## copy number 5-1000 ##
########################
prefix=tair10_cn_5_1000
min=5
max=1000

JOB=TEGS_${prefix}
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"
# Submitted batch job 10737372

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

########################
## copy number 5-2000 ##
########################
prefix=tair10_cn_5_2000
min=5
max=2000

JOB=TEGS_${prefix}
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"
# Submitted batch job 10737373

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'