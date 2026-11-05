#!/bin/bash
############
## mode 2 ##
############
# GITHUB_REPO=cloned repository of "TEgenomeSimulator"
# Navigate to the cloned repository of "TEgenomeSimulator_publication_code"
ml apptainer
ml samtools/1.20
TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # remove this line before git push
cd $WRK
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

genome=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

TIME="06:00:00"
THREADS=10
MEM=5G
cd $OUT
rm -f random_seed_tair10chr1_m2.txt
for i in {1..5};
do
prefix=tair10_ch1_m2_${i}
SEED=$SRANDOM
echo $SEED >> random_seed_tair10chr1_m2.txt
JOB=TEGS_tair10_ch1_m2_${i}
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -s $SEED -o $OUT -t $THREADS"
done
# Submitted batch job 10737314
# Submitted batch job 10737315
# Submitted batch job 10737316
# Submitted batch job 10737317
# Submitted batch job 10737318

# indexing
for i in {1..5};
do
prefix=tair10_ch1_m2_${i}
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done

##############
## mode 2+0 ##
##############
# Use user-specified TE table
ml conda
#WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # remove this line before git push
cd $WRK
#TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output

# Create chr index file
ml samtools
conda activate hraaxt_emboss # use emboss' tool infoseq to extract gc% (Version: EMBOSS:6.6.0.0) #remove hraaxt before publishing
for i in {1..5};
do
cd $OUT/TEgenomeSimulator_tair10_ch1_m2_${i}_result 
samtools faidx GCF_000001735.4_TAIR10.1_genomic_chr1.fna.masked.reformatted.nonTE.emptfixed
nontesize=$(awk '{print $2}' GCF_000001735.4_TAIR10.1_genomic_chr1.fna.masked.reformatted.nonTE.emptfixed.fai)
gc=$(cat GCF_000001735.4_TAIR10.1_genomic_chr1.fna.masked.reformatted.nonTE.emptfixed | infoseq -auto -only -name -length -pgc stdin | awk 'NR>1{print $3}')
echo "chr1,$nontesize,$gc" > $IN/tair10_chr1_index_nonte_rep${i}.csv
done
conda deactivate

cd $OUT
TIME="06:00:00"
THREADS=1
MEM=15G
rm -f random_seed_tair10chr1_m0.txt
for i in {1..5};
do
prefix=tair10_ch1_m0_${i}
SEED=$SRANDOM
echo $SEED >> random_seed_tair10chr1_m0.txt
indcsv=$IN/tair10_chr1_index_nonte_rep${i}.csv
repeat=$OUT/TEgenomeSimulator_tair10_ch1_m2_${i}_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
tetable=$OUT/TEgenomeSimulator_tair10_ch1_m2_${i}_result/TElib_sim_list_mode2.table
JOB=TEGS_tair10_ch1_m0_${1}
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 0 -p $prefix -c $indcsv -r $repeat -s $SEED -o $OUT --te_table $tetable --frag_mode 2"     
done
# Submitted batch job 10737350
# Submitted batch job 10737351
# Submitted batch job 10737352
# Submitted batch job 10737353
# Submitted batch job 10737354

# indexing
for i in {1..5};
do
prefix=tair10_ch1_m0_${i}
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done
