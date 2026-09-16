#!/bin/bash
# GITHUB_REPO=cloned repository of "TEgenomeSimulator"
# Navigate to the cloned repository of "TEgenomeSimulator_publication_code"
#TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
#WRK=/workspace/cflthc/script/TEgenomeSimulator_publication_code # remove this line before git push
cd $WRK
WRK=$(pwd)
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output/denovoTEeval
WRK=$(pwd)
mkdir -p $OUT

# copy and rename repeat fasta file
cp $GITHUB_REPO/test/input/combined_curated_TE_lib_ATOSZM_selected.fasta.gz ./
gzip -d combined_curated_TE_lib_ATOSZM_selected.fasta.gz
#repeat=/workspace/cflthc/script/TEgenomeSimulator_pub/output/tair10/output/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched # remove this line before git push
repeat=$IN/combined_curated_TE_lib_ATOSZM_selected.fasta

# prepare repeat table #
cd $IN
#cp ../output/mode0/output/TEgenomeSimulator_m0_chr10000k_result/TElib_sim_list.table ./TElib_sim_list_10000k.table
for i in {1..5}
do
cp $WRK/output/mode0/TEgenomeSimulator_m0_chr10000k_rep${i}_result/TElib_sim_list.table ./TElib_sim_list_10000k_rep${i}.table
done

# transform the format of TElib_sim_list.table to fit with the format of denovoTE-eval
for i in {1..5}
do
awk 'BEGIN{OFS="\t"}
NR==1 {
    # Print header with "#TE_id" instead of "#TE_family"
    print "#TE_id","count","idn","sd","indels","tsd","length","frag","nest"
    next
}
{
    # convert frag to integer
    frag=int($10)

    # convert tsd values
    tsd=($8=="0,0" ? "n" : "y")

    # print selected columns
    print $1,$4,$5,$6,$7,tsd,$9,frag,$11
}' TElib_sim_list_10000k_rep${i}.table > TElib_sim_list_10000k_rep${i}_denovoTEeval.table
done

# prepare config file
cd $IN
# cat ../output/mode0/output/TEgenomeSimulator_m0_chr10000k_result/TEgenomeSimulator_m0_chr10000k.yml
# prefix: "m0_chr10000k"
# rep_fasta: "../output/mode0/input/combined_curated_TE_lib_ATOSZM_selected.fasta"
# rep_list: "../output/mode0/output/TEgenomeSimulator_m0_chr10000k_result/TElib_sim_list.table"
# seed: 1
# chrs:
#   chr1:
#     prefix: "chr1"
#     seq_length: 10000000
#     gc_content: 35
for i in {1..5}
do
rm -f $OUT/random_seed_denovoTEeval.txt
SEED=$SRANDOM
echo $SEED >> $OUT/random_seed_denovoTEeval.txt
cat > config_rep${i}.yml << EOF
prefix: "m0_chr10000k_rep${i}"
prefix_nest: "m0_chr10000k_nest_rep${i}"
seq_length: 10000000
rep_fasta: "combined_curated_TE_lib_ATOSZM_selected.fasta"
rep_list: "TElib_sim_list_10000k_rep${i}_denovoTEeval.table"
seed: $SEED
gc_content: 35
EOF
done 

for i in {1..5}
do
cp -r /workspace/cflthc/bin/denovoTE-eval /workspace/cflthc/bin/denovoTE-eval_rep${i} # clean up this line before git push
cd /workspace/cflthc/bin/denovoTE-eval_rep${i} # clean up this line before git push
cp $IN/combined_curated_TE_lib_ATOSZM_selected.fasta .
cp $IN/TElib_sim_list_10000k_rep${i}_denovoTEeval.table .
cp $IN/config_rep${i}.yml config.yml
done

for i in {1..5}
do
cd /workspace/cflthc/bin/denovoTE-eval_rep${i} # clean up this line before git push
JOB=denovoTEeval_s1_rep${i}
TIME="01:00:00"
THREADS=1
MEM=10G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python ./random_sequence_TEs.py"
done
# Submitted batch job 10714048
# Submitted batch job 10714049
# Submitted batch job 10714050
# Submitted batch job 10714051
# Submitted batch job 10714052

for i in {1..5}
do
cd /workspace/cflthc/bin/denovoTE-eval_rep${i} # clean up this line before git push
JOB=denovoTEeval_s2_rep${i}
TIME="01:00:00"
THREADS=1
MEM=10G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python ./random_nest_TEs.py"
done
# Submitted batch job 10714054
# Submitted batch job 10714055
# Submitted batch job 10714056
# Submitted batch job 10714057
# Submitted batch job 10714058

# move results file to proper output directory and index
ml samtools
cd $OUT
for i in {1..5}
do
mkdir -p $OUT/denovoTEeval_rep${i}
cd $OUT/denovoTEeval_rep${i}
denovoTEeval=/workspace/cflthc/bin/denovoTE-eval_rep${i} # clean up this line before git push
mv $denovoTEeval/denovoTEeval_s* .
mv $denovoTEeval/m0_chr10000k_* .
mv $denovoTEeval/TElib_sim_list_10000k_rep${i}_denovoTEeval.table .
mv $denovoTEeval/combined_curated_TE_lib_ATOSZM_selected.fasta .
mv $denovoTEeval/config.yml .
samtools faidx m0_chr10000k_rep${i}_out_sequence_nest.fasta
samtools faidx m0_chr10000k_rep${i}_out_repeats.fasta
done


