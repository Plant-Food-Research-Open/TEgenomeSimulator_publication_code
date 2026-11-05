#!/bin/bash
#####################################################################
## Item 26 (R3-minor#10): sensitivity of the Fig. 5 recovery rate to the
## reciprocal-overlap threshold.
## The existing RepeatMasker output of the four Fig. 3D scenarios
## (replicates 1-5, 13_sim_tair10_sd_repeatmasker.sh) is re-evaluated with
## te_sim_eval_bedtools.sh at 50%, 70%, 80% and 90% reciprocal overlap.
## 80% is recomputed so that all thresholds come from the same run; it should
## match RM_recovery_eval/ (checked in 26.1_item26_overlap_threshold.Rmd).
## Output: TEgenomeSimulator_tair10_sd_scenario{n}_rep{i}_result/RM_recovery_eval_f{50,70,80,90}/
## One job per scenario x replicate; no RepeatMasker rerun is needed.
#####################################################################
# modules are loaded here and exported to every job by sbatch
ml bedtools/2.30.0
user=cflthc
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # clean up this line before git push
OUT=$WRK/output
mkdir -p $OUT/log

for n in 1 2 3 4;
do
for i in {1..5};
do
prefix=tair10_sd_scenario${n}_rep${i}
res=$OUT/TEgenomeSimulator_${prefix}_result
tegff=$res/${prefix}_repeat_annotation_out_final.gff
rmout=$res/RepeatMasker/sim_genome.fa.out
JOB=RMthreshold_${prefix}
sbatch --cpus-per-task=1 --mem 10G -t 06:00:00 -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "set -e; for f in 50 70 80 90; do mkdir -p $res/RM_recovery_eval_f\$f && cd $res/RM_recovery_eval_f\$f && bash $WRK/script/te_sim_eval_bedtools.sh $rmout $tegff 0.\$f && rm -f repeatmasker.bed truth.bed; done"
done
done
