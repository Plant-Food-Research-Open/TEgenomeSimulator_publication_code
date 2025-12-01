#!/bin/bash
############
## mode 2 ##
############
# Navigate to the cloned repository of TEgenomeSimulator
ml apptainer
WRK=$(pwd)
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT=$WRK/pub/output
mkdir -p $IN $OUT

prefix=tair10_ch1_m2
genome=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEGS_tair10_ch1_m2
TIME="06:00:00"
THREADS=10
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"

##############
## mode 2+0 ##
##############
# Navigate to the cloned repository of TEgenomeSimulator
WRK=$(pwd)
# Use user-specified TE table
ml conda
#TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT=$WRK/pub/output
prefix=tair10_ch1_m0
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta
tetable=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/TElib_sim_list_mode2.table
# Create chr index file
ml samtools
conda activate emboss # use emboss' tool infoseq to extract gc%
cd $OUT/TEgenomeSimulator_tair10_ch1_m2_result 
samtools faidx GCF_000001735.4_TAIR10.1_genomic_chr1.fna.masked.reformatted.nonTE.emptfixed
nontesize=$(awk '{print $2}' GCF_000001735.4_TAIR10.1_genomic_chr1.fna.masked.reformatted.nonTE.emptfixed.fai)
gc=$(cat GCF_000001735.4_TAIR10.1_genomic_chr1.fna.masked.reformatted.nonTE.emptfixed | infoseq -auto -only -name -length -pgc stdin | awk 'NR>1{print $3}')
echo "chr1,$nontesize,$gc" > $IN/tair10_chr1_index_nonte.csv
conda deactivate

cd $OUT
indcsv=$IN/tair10_chr1_index_nonte.csv
repeat=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
tetable=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/TElib_sim_list_mode2.table
JOB=TEGS_tair10_ch1_m0
TIME="06:00:00"
THREADS=1
MEM=15G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 0 -p $prefix -c $indcsv -r $repeat -o $OUT --te_table $tetable --frag_mode 2"     
