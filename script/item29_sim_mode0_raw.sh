#!/bin/bash

######################
## copy number 1-10 ##
######################
# GITHUB_REPO=cloned repository of "TEgenomeSimulator"
# Navigate to the cloned repository of "TEgenomeSimulator_publication_code"
ml singularity/3
ml samtools/1.20
#TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
#WRK=/workspace/cflthc/script/TEgenomeSimulator_publication_code # remove this line before git push
cd $WRK
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output/mode0
mkdir -p $IN $OUT

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
#repeat=/workspace/cflthc/script/TEgenomeSimulator_pub/output/tair10/output/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched # remove this line before git push
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
#NODE="aklppg31"  # remove this line before git push
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 0 -c $indcsv -p $prefix -r $repeat -m $max -n $min -s $SEED -o $OUT -t $THREADS"
done
done
# Submitted batch job 10704590
# Submitted batch job 10704591
# Submitted batch job 10704592
# Submitted batch job 10704593
# Submitted batch job 10704594
# Submitted batch job 10704595
# Submitted batch job 10704596
# Submitted batch job 10704597
# Submitted batch job 10704598
# Submitted batch job 10704599
# Submitted batch job 10704600
# Submitted batch job 10704601
# Submitted batch job 10704602
# Submitted batch job 10704603
# Submitted batch job 10704604
# Submitted batch job 10704605
# Submitted batch job 10704606
# Submitted batch job 10704607
# Submitted batch job 10704608
# Submitted batch job 10704609
# Submitted batch job 10704610
# Submitted batch job 10704611
# Submitted batch job 10704612
# Submitted batch job 10704613
# Submitted batch job 10704614
 
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
#NODE="aklppg31"  # remove this line before git push
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 0 -c $indcsv -p $prefix -r $repeat -m $max -n $min -s $SEED -o $OUT -t $THREADS"
done
# Submitted batch job 10704615
# Submitted batch job 10704616
# Submitted batch job 10704617
# Submitted batch job 10704618
# Submitted batch job 10704619

# indexing
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

