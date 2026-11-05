#!/bin/bash
ml apptainer/1.1

WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # clean up this line before git push
TOOL=/workspace/${user}/bin/TETools # clean up this line before git push
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT/log
cd $OUT

THREADS=20
RAM=15G
TIME="1-00:00:00"
GRP=1

for prefix in tair10_sd_scenario1 tair10_sd_scenario2 tair10_sd_scenario3 tair10_sd_scenario4
do
for i in {1..5};
do
prefix2=${prefix}_rep${i}
genome=$OUT/TEgenomeSimulator_${prefix2}_result/${prefix2}_genome_sequence_out_final.fasta
telib=$IN/athrep.updated.nonredun.noCenSatelli.fasta.stitched
JOB=RMask_${prefix2}
mkdir -p $OUT/TEgenomeSimulator_${prefix2}_result/RepeatMasker
cd $OUT/TEgenomeSimulator_${prefix2}_result/RepeatMasker
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa"
done
done
# Submitted batch job 11114090
# Submitted batch job 11114091
# Submitted batch job 11114092
# Submitted batch job 11114093
# Submitted batch job 11114094
# Submitted batch job 11114095
# Submitted batch job 11114096
# Submitted batch job 11114097
# Submitted batch job 11114098
# Submitted batch job 11114099
# Submitted batch job 11114100
# Submitted batch job 11114101
# Submitted batch job 11114102
# Submitted batch job 11114103
# Submitted batch job 11114104
# Submitted batch job 11114105
# Submitted batch job 11114106
# Submitted batch job 11114107
# Submitted batch job 11114108
# Submitted batch job 11114109

#########################################
# RepeatMasker recovery rate evaluation #
#########################################
ml bedtools/2.30.0
THREADS=1
RAM=10G
TIME="06:00:00"
for prefix in tair10_sd_scenario1 tair10_sd_scenario2 tair10_sd_scenario3 tair10_sd_scenario4
do
for i in {1..5};
do
prefix2=${prefix}_rep${i}
tegff=$OUT/TEgenomeSimulator_${prefix2}_result/${prefix2}_repeat_annotation_out_final.gff
rmout=$OUT/TEgenomeSimulator_${prefix2}_result/RepeatMasker/sim_genome.fa.out
# bedtools intersect
out=$OUT/TEgenomeSimulator_${prefix2}_result/RM_recovery_eval
mkdir -p $out
cd $out
JOB=RMCOVeval_${prefix2}
cut=0.8
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$WRK/script/te_sim_eval_bedtools.sh $rmout $tegff $cut"
done
done
# Submitted batch job 11114126
# Submitted batch job 11114127
# Submitted batch job 11114128
# Submitted batch job 11114129
# Submitted batch job 11114130
# Submitted batch job 11114131
# Submitted batch job 11114132
# Submitted batch job 11114133
# Submitted batch job 11114134
# Submitted batch job 11114135
# Submitted batch job 11114136
# Submitted batch job 11114137
# Submitted batch job 11114138
# Submitted batch job 11114139
# Submitted batch job 11114140
# Submitted batch job 11114141
# Submitted batch job 11114142
# Submitted batch job 11114143
# Submitted batch job 11114144
# Submitted batch job 11114145



