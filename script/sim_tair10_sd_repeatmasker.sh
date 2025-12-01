#!/bin/bash
ml singularity/3
WRK=$(pwd)
TOOL=$(cd ../../bin/TETools && pwd)
IN=$WRK/pub/input
OUT=../TEgenomeSimulator_pub/output
cd $OUT
OUT=$(pwd)

THREADS=20
RAM=15G
TIME="1-00:00:00"
GRP=1

for prefix in tair10_sd_scenario1 tair10_sd_scenario2 tair10_sd_scenario3 tair10_sd_scenario4
do
genome=$OUT/tair10/output/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
telib=$OUT/tair10/output/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
JOB=RMask_$prefix
mkdir -p $OUT/tair10/output/TEgenomeSimulator_${prefix}_result/RepeatMasker
cd $OUT/tair10/output/TEgenomeSimulator_${prefix}_result/RepeatMasker
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa"
done

#########################################
# RepeatMasker recovery rate evaluation #
#########################################
ml bedtools/2.30.0
THREADS=1
RAM=10G
TIME="06:00:00"
for prefix in tair10_sd_scenario1 tair10_sd_scenario2 tair10_sd_scenario3 tair10_sd_scenario4
do
tegff=$OUT/tair10/output/TEgenomeSimulator_${prefix}_result/${prefix}_repeat_annotation_out_final.gff
rmout=$OUT/tair10/output/TEgenomeSimulator_${prefix}_result/RepeatMasker/sim_genome.fa.out
# bedtools intersect
out=$OUT/tair10/output/TEgenomeSimulator_${prefix}_result/RM_recovery_eval
mkdir -p $out
cd $out
JOB=RMCOVeval_${prefix}
cut=0.8
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$WRK/pub/script/te_sim_eval_bedtools.sh $rmout $tegff $cut"
done





