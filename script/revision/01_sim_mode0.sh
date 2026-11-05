#!/bin/bash

##################################
## TEGS mode 0 copy number 1-10 ##
##################################
# GITHUB_REPO=cloned repository of "TEgenomeSimulator"
# Navigate to the cloned repository of "TEgenomeSimulator_publication_code"
ml apptainer
ml samtools/1.20
TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # update this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # clean up this line before git push
cd $WRK
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output/mode0
mkdir -p $IN $OUT/log

cd $IN
echo "chr1,10000,35" > random_genome_chr_index_10k.csv
echo "chr1,100000,35" > random_genome_chr_index_100k.csv
echo "chr1,1000000,35" > random_genome_chr_index_1000k.csv
echo "chr1,10000000,35" > random_genome_chr_index_10000k.csv
echo "chr1,100000000,35" > random_genome_chr_index_100000k.csv

echo "chr1,100000000,35" > random_genome_chr_index_1000000k.csv
echo "chr2,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr3,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr4,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr5,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr6,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr7,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr8,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr9,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr10,100000000,35" >> random_genome_chr_index_1000000k.csv

cp $GITHUB_REPO/test/input/combined_curated_TE_lib_ATOSZM_selected.fasta.gz ./
gzip -d combined_curated_TE_lib_ATOSZM_selected.fasta.gz
repeat=$IN/combined_curated_TE_lib_ATOSZM_selected.fasta

min=1
max=10

for chr in 10 100 1000 10000 100000
do 
rm -f random_seed_m0_chr${chr}k.txt
indcsv=$IN/random_genome_chr_index_${chr}k.csv
for i in {1..5};
do
SEED=$SRANDOM
echo $SEED >> random_seed_m0_chr${chr}k.txt
prefix=m0_chr${chr}k_rep${i}
JOB=TEgenomeSimulator_m0_chr${chr}k_${i}
TIME="1-00:00:00"
THREADS=1
MEM=30G
NODE="aklppg31"  # remove this line before git push
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 0 -c $indcsv -p $prefix -r $repeat -m $max -n $min -s $SEED -o $OUT -t $THREADS"
done
done
# Submitted batch job 10733351
# Submitted batch job 10733352
# Submitted batch job 10733353
# Submitted batch job 10733354
# Submitted batch job 10733355
# Submitted batch job 10733356
# Submitted batch job 10733357
# Submitted batch job 10733358
# Submitted batch job 10733359
# Submitted batch job 10733360
# Submitted batch job 10733361
# Submitted batch job 10733362
# Submitted batch job 10733363
# Submitted batch job 10733364
# Submitted batch job 10733365
# Submitted batch job 10733366
# Submitted batch job 10733367
# Submitted batch job 10733368
# Submitted batch job 10733369
# Submitted batch job 10733370
# Submitted batch job 10733371
# Submitted batch job 10733372
# Submitted batch job 10733373
# Submitted batch job 10733374
# Submitted batch job 10733375

# 20261012 for FigS1 (chr=10000)
# Submitted batch job 10956520
# Submitted batch job 10956521
# Submitted batch job 10956522
# Submitted batch job 10956523
# Submitted batch job 10956524

chr=1000000
indcsv=$IN/random_genome_chr_index_${chr}k.csv
rm -f random_seed_m0_chr${chr}k.txt
for i in {1..5};
do
SEED=$SRANDOM
echo $SEED >> random_seed_m0_chr${chr}k.txt
prefix=m0_chr${chr}k_rep${i}
JOB=TEgenomeSimulator_m0_chr${chr}k_${i}
TIME="1-00:00:00"
THREADS=1
MEM=100G
NODE="aklppg31"  # remove this line before git push
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 0 -c $indcsv -p $prefix -r $repeat -m $max -n $min -s $SEED -o $OUT -t $THREADS"
done
# Submitted batch job 10733378
# Submitted batch job 10733379
# Submitted batch job 10733380
# Submitted batch job 10733381
# Submitted batch job 10733382

##################################
## indexing
##################################
for chr in 10 100 1000 10000 100000
do 
for i in {1..5};
do
prefix=m0_chr${chr}k_rep${i}
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done
done

chr=1000000
for i in {1..5};
do
prefix=m0_chr${chr}k_rep${i}
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done

##################################
## RepeatMasker for 10MB sim genome (10000k)
## Using TEgenomeSimulator's mode 2 with --te_scan_only
##################################
ml singularity/3
WRK=$(pwd)
#TOOL=$(cd ../../bin/TETools && pwd)
IN=$WRK/input
OUT=$WRK/output/mode0
cd $OUT

THREADS=20
MEM=15G
TIME="1-00:00:00"
for i in {1..5};
do
prefix=m0_chr10000k_rep${i}
genome=$OUT/TEgenomeSimulator_$prefix'_result'/m0_chr10000k_rep${i}_genome_sequence_out_final.fasta
repeat=$IN/combined_curated_TE_lib_ATOSZM_selected.fasta
JOB=te_scan_only_$prefix
mkdir -p $OUT/TEgenomeSimulator_$prefix'_result'/te_scan_only
cd $OUT/TEgenomeSimulator_$prefix'_result'/te_scan_only
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT/TEgenomeSimulator_$prefix'_result'/te_scan_only -t $THREADS --te_scan_only"
done
# Submitted batch job 10956525
# Submitted batch job 10956526
# Submitted batch job 10956527
# Submitted batch job 10956528
# Submitted batch job 10956529