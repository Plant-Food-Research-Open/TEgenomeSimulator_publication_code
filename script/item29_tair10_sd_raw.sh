#!/bin/bash
##########################
## scenario I 
## high idn 90-100 
## high sd 15-20
##########################
ml singularity/3
ml samtools/1.20
#TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
#WRK=/workspace/cflthc/script/TEgenomeSimulator_publication_code # remove this line before git push
#TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
cd $WRK
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

minidn=90
maxidn=100
minsd=15
maxsd=20
min=5
max=1000
#genome=/workspace/cflthc/script/TEgenomeSimulator_pub/output/tair10/output/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed # remove this line before git push
#repeat=/workspace/cflthc/script/TEgenomeSimulator_pub/output/tair10/output/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched # remove this line before git push
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

TIME="1-00:00:00"
THREADS=1
MEM=10G
#NODE="aklppg31"  # remove this line before git push
cd $OUT
rm -f random_seed_tair10_sd_scenario1.txt
for i in {1..5};
do
SEED=$SRANDOM
echo $SEED >> random_seed_tair10_sd_scenario1.txt
JOB=TEgenomeSimulator_scenario1_${i}
prefix=tair10_sd_scenario1_rep${i}
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -s $SEED -o $OUT"
done
# Submitted batch job 10704521
# Submitted batch job 10704522
# Submitted batch job 10704523
# Submitted batch job 10704524
# Submitted batch job 10704525

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
ml singularity/3
ml samtools/1.20
#TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
#WRK=/workspace/cflthc/script/TEgenomeSimulator_publication_code # remove this line before git push
#TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

minidn=70
maxidn=80
minsd=15
maxsd=20
min=5
max=1000
#genome=/workspace/cflthc/script/TEgenomeSimulator_pub/output/tair10/output/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed # remove this line before git push
#repeat=/workspace/cflthc/script/TEgenomeSimulator_pub/output/tair10/output/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched # remove this line before git push
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

TIME="1-00:00:00"
THREADS=1
MEM=10G
#NODE="aklppg31"  # remove this line before git push
cd $OUT
rm -f random_seed_tair10_sd_scenario2.txt
for i in {1..5};
do
SEED=$SRANDOM
echo $SEED >> random_seed_tair10_sd_scenario2.txt
JOB=TEgenomeSimulator_scenario2_${i}
prefix=tair10_sd_scenario2_rep${i}
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -s $SEED -o $OUT"
done
# Submitted batch job 10704526
# Submitted batch job 10704527
# Submitted batch job 10704528
# Submitted batch job 10704529
# Submitted batch job 10704530

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
ml singularity/3
ml samtools/1.20
#TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
#WRK=/workspace/cflthc/script/TEgenomeSimulator_publication_code # remove this line before git push
#TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

minidn=70
maxidn=80
minsd=1
maxsd=5
min=5
max=1000
#genome=/workspace/cflthc/script/TEgenomeSimulator_pub/output/tair10/output/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed # remove this line before git push
#repeat=/workspace/cflthc/script/TEgenomeSimulator_pub/output/tair10/output/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched # remove this line before git push
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

TIME="1-00:00:00"
THREADS=1
MEM=10G
#NODE="aklppg31"  # remove this line before git push
cd $OUT
rm -f random_seed_tair10_sd_scenario3.txt
for i in {1..5};
do
SEED=$SRANDOM
echo $SEED >> random_seed_tair10_sd_scenario3.txt
JOB=TEgenomeSimulator_scenario3_${i}
prefix=tair10_sd_scenario3_rep${i}
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -s $SEED -o $OUT"
done
# Submitted batch job 10704532
# Submitted batch job 10704533
# Submitted batch job 10704534
# Submitted batch job 10704535
# Submitted batch job 10704536

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
ml singularity/3
ml samtools/1.20
#TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
#WRK=/workspace/cflthc/script/TEgenomeSimulator_publication_code # remove this line before git push
#TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

minidn=90
maxidn=100
minsd=1
maxsd=5
min=5
max=1000
#genome=/workspace/cflthc/script/TEgenomeSimulator_pub/output/tair10/output/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed # remove this line before git push
#repeat=/workspace/cflthc/script/TEgenomeSimulator_pub/output/tair10/output/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched # remove this line before git push
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

TIME="1-00:00:00"
THREADS=1
MEM=10G
#NODE="aklppg31"  # remove this line before git push
cd $OUT
rm -f random_seed_tair10_sd_scenario4.txt
for i in {1..5};
do
SEED=$SRANDOM
echo $SEED >> random_seed_tair10_sd_scenario4.txt
JOB=TEgenomeSimulator_scenario4_${i}
prefix=tair10_sd_scenario4_rep${i}
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -s $SEED -o $OUT"
done
# Submitted batch job 10704537
# Submitted batch job 10704538
# Submitted batch job 10704539
# Submitted batch job 10704540
# Submitted batch job 10704541

# indexing
for i in {1..5};
do
prefix=tair10_sd_scenario4_rep${i}
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done