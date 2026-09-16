#!/bin/bash
# Navigate to the cloned repository of "TEgenomeSimulator_publication_code"

ml conda
#conda create --name tesimeval python=3.9 
conda activate cflthc_tesimeval # clean up this line before git push
#pip install biopython numpy scipy scikit-learn pandas matplotlib
#WRK=/workspace/cflthc/script/TEgenomeSimulator_publication_code # remove this line before git push
cd $WRK
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output

################################################
# Run te_sim_eval.py 
################################################

#### TEGSm2 ####
for i in {1..5}
do
prefix=tair10_ch1_m2_${i}
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
tegff=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_repeat_annotation_out_final.gff
THREADS=1
RAM=30G
TIME="6-00:00:00"
for eval in kmer entropy
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/script/te_sim_eval.py --genome $genome --te_gff $tegff --outdir $out --k 6 --win 1000 --eval ${eval}"
done
done
# Submitted batch job 10722091
# Submitted batch job 10722092
# Submitted batch job 10722093
# Submitted batch job 10722094
# Submitted batch job 10722095
# Submitted batch job 10722096
# Submitted batch job 10722097
# Submitted batch job 10722098
# Submitted batch job 10722099
# Submitted batch job 10722100

#### TEGSm2+0 ####
for i in {1..5}
do
prefix=tair10_ch1_m0_${i}
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
tegff=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_repeat_annotation_out_final.gff
THREADS=1
RAM=30G
TIME="6-00:00:00"
for eval in kmer entropy
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/script/te_sim_eval.py --genome $genome --te_gff $tegff --outdir $out --k 6 --win 1000 --eval ${eval}"
done
done
# Submitted batch job 10722121
# Submitted batch job 10722122
# Submitted batch job 10722123
# Submitted batch job 10722124
# Submitted batch job 10722125
# Submitted batch job 10722126
# Submitted batch job 10722127
# Submitted batch job 10722128
# Submitted batch job 10722129
# Submitted batch job 10722130
#### Garlic ####
for i in {1..5}
do
prefix=Garlic_rep${i}
genome=$OUT/Garlic/${prefix}/tair10_garlic_sim_chr1.fa.fasta
tegff=$OUT/Garlic/${prefix}/tair10_garlic_sim_chr1.fa.inserts.te.gff
THREADS=1
RAM=30G
TIME="6-00:00:00"
for eval in kmer entropy
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/script/te_sim_eval.py --genome $genome --te_gff $tegff --outdir $out --k 6 --win 1000 --eval ${eval}"
done
done
# Submitted batch job 10722135
# Submitted batch job 10722136
# Submitted batch job 10722137
# Submitted batch job 10722138
# Submitted batch job 10722139
# Submitted batch job 10722140
# Submitted batch job 10722141
# Submitted batch job 10722142
# Submitted batch job 10722143
# Submitted batch job 10722144

#### Original #### 
prefix=original
genome=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
tegff=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out.gff

THREADS=1
RAM=30G
TIME="6-00:00:00"

for eval in kmer entropy 
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/script/te_sim_eval.py --genome $genome --te_gff $tegff --outdir $out --k 6 --win 1000 --eval ${eval}"
done
# Submitted batch job 10722146
# Submitted batch job 10722147

