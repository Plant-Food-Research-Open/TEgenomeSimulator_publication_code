#!/bin/bash
# GITHUB_REPO=cloned repository of "TEgenomeSimulator"
# Navigate to the cloned repository of "TEgenomeSimulator_publication_code"
TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # clean up this line before git push
cd $WRK
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output/denovoTEeval
mkdir -p $OUT

# copy and rename repeat fasta file
cp $GITHUB_REPO/test/input/combined_curated_TE_lib_ATOSZM_selected.fasta.gz ./
gzip -d combined_curated_TE_lib_ATOSZM_selected.fasta.gz
repeat=$IN/combined_curated_TE_lib_ATOSZM_selected.fasta

# prepare repeat table
cd $IN
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

TOOL_DIR=/workspace/${user}/bin # The directory where denovoTE-eval cloned into
for i in {1..5}
do
cp -r $TOOL_DIR/denovoTE-eval $TOOL_DIR/denovoTE-eval_rep${i}
cd $TOOL_DIR/denovoTE-eval_rep${i} 
cp $IN/combined_curated_TE_lib_ATOSZM_selected.fasta .
cp $IN/TElib_sim_list_10000k_rep${i}_denovoTEeval.table .
cp $IN/config_rep${i}.yml config.yml
done

for i in {1..5}
do
cd $TOOL_DIR/denovoTE-eval_rep${i} 
JOB=denovoTEeval_s1_rep${i}
TIME="01:00:00"
THREADS=1
MEM=10G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python ./random_sequence_TEs.py"
done
# Submitted batch job 10737184
# Submitted batch job 10737185
# Submitted batch job 10737186
# Submitted batch job 10737187
# Submitted batch job 10737188

for i in {1..5}
do
cd $TOOL_DIR/denovoTE-eval_rep${i}
JOB=denovoTEeval_s2_rep${i}
TIME="01:00:00"
THREADS=1
MEM=10G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python ./random_nest_TEs.py"
done
# Submitted batch job 10737189
# Submitted batch job 10737190
# Submitted batch job 10737191
# Submitted batch job 10737192
# Submitted batch job 10737193

# move results file to proper output directory and index
ml samtools
cd $OUT
for i in {1..5}
do
mkdir -p $OUT/denovoTEeval_rep${i}
cd $OUT/denovoTEeval_rep${i}
denovoTEeval=$TOOL_DIR/denovoTE-eval_rep${i}
mv $denovoTEeval/denovoTEeval_s* .
mv $denovoTEeval/m0_chr10000k_* .
mv $denovoTEeval/TElib_sim_list_10000k_rep${i}_denovoTEeval.table .
mv $denovoTEeval/combined_curated_TE_lib_ATOSZM_selected.fasta .
mv $denovoTEeval/config.yml .
samtools faidx m0_chr10000k_rep${i}_out_sequence_nest.fasta
samtools faidx m0_chr10000k_rep${i}_out_repeats.fasta
done


##################################
## RepeatMasker for 10MB sim genome (10000k)
## Using TEgenomeSimulator's mode 2 with --te_scan_only
##################################
ml singularity/3
WRK=$(pwd)
#TOOL=$(cd ../../bin/TETools && pwd)
IN=$WRK/input
OUT=$WRK/output/denovoTEeval
#OUT=$WRK/output/mode0
cd $OUT

THREADS=20
MEM=15G
TIME="1-00:00:00"
for i in {1..5};
do
prefix=denovoTEeval_rep${i}
genome=$OUT/${prefix}/m0_chr10000k_rep${i}_out_sequence_nest.fasta
repeat=$IN/combined_curated_TE_lib_ATOSZM_selected.fasta
JOB=te_scan_only_$prefix
mkdir -p $OUT/${prefix}/te_scan_only
cd $OUT/${prefix}/te_scan_only
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT/${prefix}/te_scan_only -t $THREADS --te_scan_only"
done
# Submitted batch job 10956540
# Submitted batch job 10956541
# Submitted batch job 10956542
# Submitted batch job 10956543
# Submitted batch job 10956544