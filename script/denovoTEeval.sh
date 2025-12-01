#!/bin/bash
WRK=$(pwd)
IN=$WRK/pub/input
OUT=$WRK/pub/output/denovoTEeval
mkdir -p $IN $OUT

# copy and rename repeat fasta file
cp ../../test/input/combined_curated_TE_lib_ATOSZM_selected.fasta.gz .
gzip -d -k combined_curated_TE_lib_ATOSZM_selected.fasta.gz

# prepare repeat table #
cd $IN
cp ../output/mode0/output/TEgenomeSimulator_m0_chr10000k_result/TElib_sim_list.table ./TElib_sim_list_10000k.table

# transform the format of TElib_sim_list.table to fit with the format of denovoTE-eval
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
}' TElib_sim_list_10000k.table > TElib_sim_list_10000k_denovoTEeval.table

# prepare config file
cd $IN
cat ../output/mode0/output/TEgenomeSimulator_m0_chr10000k_result/TEgenomeSimulator_m0_chr10000k.yml
# prefix: "m0_chr10000k"
# rep_fasta: "../output/mode0/input/combined_curated_TE_lib_ATOSZM_selected.fasta"
# rep_list: "../output/mode0/output/TEgenomeSimulator_m0_chr10000k_result/TElib_sim_list.table"
# seed: 1
# chrs:
#   chr1:
#     prefix: "chr1"
#     seq_length: 10000000
#     gc_content: 35
cat > config.yml << EOF
prefix: "m0_chr10000k"
prefix_nest: "m0_chr10000k_nest"
seq_length: 10000000
rep_fasta: "combined_curated_TE_lib_ATOSZM_selected.fasta"
rep_list: "TElib_sim_list_10000k_denovoTEeval.table"
seed: 1
gc_content: 35
EOF

cd ../../../../bin/denovoTE-eval
cp ../../script/TEgenomeSimulator/pub/input/combined_curated_TE_lib_ATOSZM_selected.fasta .
cp ../../script/TEgenomeSimulator/pub/input/TElib_sim_list_10000k_denovoTEeval.table .
cp ../../script/TEgenomeSimulator/pub/input/config.yml .

JOB=denovoTEeval_s1
TIME="01:00:00"
THREADS=1
MEM=10G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python ./random_sequence_TEs.py"

JOB=denovoTEeval_s2
TIME="01:00:00"
THREADS=1
MEM=10G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python ./random_nest_TEs.py"

# move results file to proper output directory
cd $OUT
denovoTEeval='../../../../../bin/denovoTE-eval'
mv $denovoTEeval/denovoTEeval_s* .
mv $denovoTEeval/m0_chr10000k_out_* .
mv $denovoTEeval/TElib_sim_list_denovoTEeval.table .
mv $denovoTEeval/TElib_sim_list_10000k_denovoTEeval.table .
mv $denovoTEeval/combined_curated_TE_lib_ATOSZM_selected.fasta .
mv $denovoTEeval/config.yml .

ml samtools
samtools faidx m0_chr10000k_out_sequence_nest.fasta
samtools faidx m0_chr10000k_out_repeats.fasta