#!/bin/bash
#####################################################################
## Items 34 & 35: sensitivity of RepeatMasker recovery to the nested
## insertion rate and to the beta (integrity) distribution
##
## Base: Fig. 3D scenario 1 (identity 90-100, SD 15-20), replicates 1-3.
## Each new run reuses that replicate's TE mutagenesis table
## (TElib_sim_list.table) and seed, so copy numbers, identities and the
## random stream match the baseline; only one setting changes per config:
##   nest0       - nest column = 0 for every family
##   nest30      - nest column = 30 for every LTR/Copia and LTR/Gypsy family
##   fragheavy   - beta(alpha = 0.5, beta = 2)   (mostly fragmented)
##   intactheavy - beta(alpha = 2,   beta = 0.5) (mostly intact)
## The baseline (random nest 0-30 per family, beta(0.5, 0.7)) is the
## existing tair10_sd_scenario1_rep1-3 runs, which are not rerun.
##
## Each config/replicate is submitted as a chain of three jobs:
## simulation -> RepeatMasker -> recovery evaluation (afterok dependencies)
#####################################################################
# modules are loaded here and exported to every job by sbatch
ml apptainer/1.1
ml samtools/1.20
ml bedtools/2.30.0
user=cflthc
TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # update this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # clean up this line before git push
TOOL=/workspace/${user}/bin/TETools # clean up this line before git push
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT/log

genome=$OUT/TEgenomeSimulator_tair10_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
telib=$IN/athrep.updated.nonredun.noCenSatelli.fasta.stitched
TABLES=$OUT/tair10_sensitivity_tables
mkdir -p $TABLES

##########################
## 1. TE tables per replicate and config
##########################
for i in {1..3};
do
base=$OUT/TEgenomeSimulator_tair10_sd_scenario1_rep${i}_result/TElib_sim_list.table
# baseline table, used unchanged by fragheavy/intactheavy
cp $base $TABLES/scenario1_rep${i}_base.table
# column 3 = superfamily, column 11 = nest (% of the family's copies inserted as nested)
awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} {$11=0; print}' $base > $TABLES/scenario1_rep${i}_nest0.table
awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} {$11=(toupper($3) ~ /^LTR\/(COPIA|GYPSY)/ ? 30 : 0); print}' $base > $TABLES/scenario1_rep${i}_nest30.table
done

##########################
## 2. Simulation -> RepeatMasker -> recovery evaluation
##########################
# config name, TE table suffix, alpha, beta
configs=(
"nest0 nest0 0.5 0.7"
"nest30 nest30 0.5 0.7"
"fragheavy base 0.5 2"
"intactheavy base 2 0.5"
)

for i in {1..3};
do
# same seed as the baseline replicate
SEED=$(awk '$1=="seed:"{print $2}' $OUT/TEgenomeSimulator_tair10_sd_scenario1_rep${i}_result/TEgenomeSimulator_tair10_sd_scenario1_rep${i}.yml)
for cfg in "${configs[@]}";
do
read name table_suffix alpha beta <<< "$cfg"
prefix=tair10_sens_${name}_rep${i}
res=$OUT/TEgenomeSimulator_${prefix}_result
te_table=$TABLES/scenario1_rep${i}_${table_suffix}.table

# simulation (+ faidx)
JOB=TEgenomeSimulator_sens_${name}_${i}
jid_sim=$(sbatch --parsable --cpus-per-task=1 --mem 10G -t 1-00:00:00 -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat --te_table $te_table -a $alpha -b $beta -s $SEED -o $OUT && cd $res && samtools faidx ${prefix}_genome_sequence_out_final.fasta && samtools faidx ${prefix}_repeat_sequence_out_final.fasta")

# RepeatMasker (same settings as 13_sim_tair10_sd_repeatmasker.sh)
mkdir -p $res/RepeatMasker
JOB=RMask_${prefix}
jid_rm=$(sbatch --parsable --dependency=afterok:$jid_sim --cpus-per-task=20 --mem 15G -t 1-00:00:00 -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "cd $res/RepeatMasker && ln -sf $res/${prefix}_genome_sequence_out_final.fasta sim_genome.fa && ln -sf $telib sim_used_telib.fa && $TOOL/dfam-tetools.sif RepeatMasker -pa 20 -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa")

# recovery evaluation (80% reciprocal overlap)
mkdir -p $res/RM_recovery_eval
JOB=RMCOVeval_${prefix}
jid_eval=$(sbatch --parsable --dependency=afterok:$jid_rm --cpus-per-task=1 --mem 10G -t 06:00:00 -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "cd $res/RM_recovery_eval && $WRK/script/te_sim_eval_bedtools.sh $res/RepeatMasker/sim_genome.fa.out $res/${prefix}_repeat_annotation_out_final.gff 0.8")

echo "$prefix  seed=$SEED  sim=$jid_sim  rm=$jid_rm  eval=$jid_eval"
done
done
