#!/bin/bash
##########################
## scenario I 
## high idn 90-100 
## high sd 15-20
##########################
ml apptainer/1.1
ml samtools/1.20

TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # update this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # clean up this line before git push
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT/log

minidn=90
maxidn=100
minsd=15
maxsd=20
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

TIME="1-00:00:00"
THREADS=1
MEM=10G
cd $OUT
rm -f random_seed_tair10_sd_scenario1.txt
for i in {1..5};
do
SEED=$SRANDOM
echo $SEED >> random_seed_tair10_sd_scenario1.txt
JOB=TEgenomeSimulator_scenario1_${i}
prefix=tair10_sd_scenario1_rep${i}
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -s $SEED -o $OUT"
done
# Submitted batch job 10957565
# Submitted batch job 10957566
# Submitted batch job 10957567
# Submitted batch job 10957568
# Submitted batch job 10957569

# indexing
for i in {1..5};
do
prefix=tair10_sd_scenario1_rep${i}
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done

##########################
## scenario II
## low idn 70-80 
## high sd 15-20
##########################
ml apptainer/1.1
ml samtools/1.20
TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # update this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # clean up this line before git push
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT/log

minidn=70
maxidn=80
minsd=15
maxsd=20
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

TIME="1-00:00:00"
THREADS=1
MEM=10G
cd $OUT
rm -f random_seed_tair10_sd_scenario2.txt
for i in {1..5};
do
SEED=$SRANDOM
echo $SEED >> random_seed_tair10_sd_scenario2.txt
JOB=TEgenomeSimulator_scenario2_${i}
prefix=tair10_sd_scenario2_rep${i}
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -s $SEED -o $OUT"
done
# Submitted batch job 10957570
# Submitted batch job 10957571
# Submitted batch job 10957572
# Submitted batch job 10957573
# Submitted batch job 10957574

# indexing
for i in {1..5};
do
prefix=tair10_sd_scenario2_rep${i}
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done

##########################
## scenario III
## low idn 70-80 
## low sd 1-5
##########################
ml apptainer/1.1
ml samtools/1.20
TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # update this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # clean up this line before git push
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT/log

minidn=70
maxidn=80
minsd=1
maxsd=5
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

TIME="1-00:00:00"
THREADS=1
MEM=10G
cd $OUT
rm -f random_seed_tair10_sd_scenario3.txt
for i in {1..5};
do
SEED=$SRANDOM
echo $SEED >> random_seed_tair10_sd_scenario3.txt
JOB=TEgenomeSimulator_scenario3_${i}
prefix=tair10_sd_scenario3_rep${i}
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -s $SEED -o $OUT"
done
# Submitted batch job 10957575
# Submitted batch job 10957576
# Submitted batch job 10957577
# Submitted batch job 10957578
# Submitted batch job 10957579

# indexing
for i in {1..5};
do
prefix=tair10_sd_scenario3_rep${i}
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done

##########################
## scenario IV
## high idn 90-100 
## low sd 1-5
##########################
ml apptainer/1.1
ml samtools/1.20
TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # update this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # clean up this line before git push
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT/log

minidn=90
maxidn=100
minsd=1
maxsd=5
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

TIME="1-00:00:00"
THREADS=1
MEM=10G
cd $OUT
rm -f random_seed_tair10_sd_scenario4.txt
for i in {1..5};
do
SEED=$SRANDOM
echo $SEED >> random_seed_tair10_sd_scenario4.txt
JOB=TEgenomeSimulator_scenario4_${i}
prefix=tair10_sd_scenario4_rep${i}
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -s $SEED -o $OUT"
done
# Submitted batch job 10957580
# Submitted batch job 10957581
# Submitted batch job 10957582
# Submitted batch job 10957583
# Submitted batch job 10957584

# indexing
for i in {1..5};
do
prefix=tair10_sd_scenario4_rep${i}
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done