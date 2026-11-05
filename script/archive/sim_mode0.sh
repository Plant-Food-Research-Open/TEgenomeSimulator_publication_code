#!/bin/bash

######################
## copy number 1-10 ##
######################
# Navigate to the cloned repository of TEgenomeSimulator
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/pub/input
OUT=$WRK/pub/output/mode0
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

cp ../../test/input/combined_curated_TE_lib_ATOSZM_selected.fasta.gz ./
gzip -d combined_curated_TE_lib_ATOSZM_selected.fasta.gz
repeat=$IN/combined_curated_TE_lib_ATOSZM_selected.fasta

min=1
max=10

for chr in 10 100 1000 10000 100000
do 
indcsv=$IN/random_genome_chr_index_${chr}k.csv
prefix=m0_chr${chr}k
JOB=TEgenomeSimulator_m0_chr${chr}k
TIME="1-00:00:00"
THREADS=1
MEM=30G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 0 -c $indcsv -p $prefix -r $repeat -m $max -n $min -o $OUT -t $THREADS"
done

chr=1000000
indcsv=$IN/random_genome_chr_index_${chr}k.csv
prefix=m0_chr${chr}k
JOB=TEgenomeSimulator_m0_chr${chr}k
TIME="1-00:00:00"
THREADS=1
MEM=100G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 0 -c $indcsv -p $prefix -r $repeat -m $max -n $min -o $OUT -t $THREADS"

# indexing
for chr in 10 100 1000 10000 100000
do 
prefix=m0_chr${chr}k
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done

chr=1000000
prefix=m0_chr${chr}k
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
