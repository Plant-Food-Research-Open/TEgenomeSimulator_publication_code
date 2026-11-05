#!/bin/bash
#####################################################################
## Items 23-28: full-metric benchmark of RepeatMasker on scenario 1
## (Fig. 3D, identity 90-100, SD 15-20), replicates 1-5, with
##   full   - the complete library used for simulation (existing RepeatMasker runs)
##   drop25 - 25% of library families removed (random, fixed seed)
##   drop50 - 50% removed (a superset of drop25's removed families)
## Each annotation is scored against the simulator's ground-truth GFF with
## TEbenchtop mode 2 (bp-level metrics, by TE type and identity x integrity
## slice) and mode 3 (boundaries, coverage), plus false-positive loci.
## RepeatMasker is also run on the TE-depleted backbone alone, to measure
## how many false positives come from TE remnants left in the backbone.
##
## Step 1 runs immediately; step 2 submits RepeatMasker jobs; step 3
## submits one TEbenchtop head job per replicate x library, chained to the
## RepeatMasker job where needed (afterok). Run the whole script once.
## To resubmit only step 3 (RepeatMasker already run or still running):
##   STEP3_ONLY=1 bash script/revision/15_tair10_benchmark_repeatmasker.sh
#####################################################################
# modules are loaded here and exported to every job by sbatch
ml apptainer/1.1
ml bedtools/2.30.0
user=cflthc
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # clean up this line before git push
TOOL=/workspace/${user}/bin/TETools # clean up this line before git push
TEBENCHTOP=/workspace/${user}/script/TEbenchtop # clean up this line before git push
IN=$WRK/input
OUT=$WRK/output
mkdir -p $OUT/log

telib=$IN/athrep.updated.nonredun.noCenSatelli.fasta.stitched
backbone=$OUT/TEgenomeSimulator_tair10_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
LIBDIR=$OUT/tair10_benchmark_libraries
mkdir -p $LIBDIR

declare -A jid_rm
if [[ "${STEP3_ONLY:-0}" != "1" ]]; then
##########################
## 1. Reduced libraries
##########################
# one random order of families (fixed seed); drop25 removes the first 25%, drop50 the first 50%
grep '>' $telib | sed 's/^>//' |\
awk -v seed=20261104 'BEGIN{srand(seed)} {printf "%.10f\t%s\n", rand(), $0}' |\
sort -k1,1g | cut -f2 > $LIBDIR/family_random_order.txt
nfam=$(wc -l < $LIBDIR/family_random_order.txt)
for pct in 25 50;
do
head -n $(( nfam * pct / 100 )) $LIBDIR/family_random_order.txt > $LIBDIR/dropped_families_drop${pct}.txt
awk 'NR==FNR{drop[">"$0]=1; next} /^>/{keep=!($0 in drop)} keep' $LIBDIR/dropped_families_drop${pct}.txt $telib > $LIBDIR/telib_drop${pct}.fa
echo "drop${pct}: $(wc -l < $LIBDIR/dropped_families_drop${pct}.txt) of $nfam families removed; $(grep -c '>' $LIBDIR/telib_drop${pct}.fa) kept"
done

##########################
## 2. RepeatMasker (same settings as 13_sim_tair10_sd_repeatmasker.sh)
##########################
run_repeatmasker () { # $1 = run directory, $2 = genome fasta, $3 = library, $4 = job name; prints the job id
  mkdir -p $1
  sbatch --parsable --cpus-per-task=20 --mem 15G -t 1-00:00:00 -J $4 -o $OUT/log/$4.out -e $OUT/log/$4.err --wrap "cd $1 && ln -sf $2 sim_genome.fa && ln -sf $3 sim_used_telib.fa && $TOOL/dfam-tetools.sif RepeatMasker -pa 20 -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa"
}

for i in {1..5};
do
res=$OUT/TEgenomeSimulator_tair10_sd_scenario1_rep${i}_result
for lib in drop25 drop50;
do
jid_rm[${i}_${lib}]=$(run_repeatmasker $res/RepeatMasker_${lib} $res/tair10_sd_scenario1_rep${i}_genome_sequence_out_final.fasta $LIBDIR/telib_${lib}.fa RMask_scenario1_rep${i}_${lib})
echo "scenario1_rep${i} ${lib}: RepeatMasker job ${jid_rm[${i}_${lib}]}"
done
done

# backbone only, full library (false-positive baseline)
jid_backbone=$(run_repeatmasker $OUT/TEgenomeSimulator_tair10_result/RepeatMasker_backbone $backbone $telib RMask_tair10_backbone)
sbatch --dependency=afterok:$jid_backbone --cpus-per-task=1 --mem 2G -t 00:30:00 -J RMgff_backbone -o $OUT/log/RMgff_backbone.out -e $OUT/log/RMgff_backbone.err --wrap "cd $OUT/TEgenomeSimulator_tair10_result/RepeatMasker_backbone && bash $TEBENCHTOP/bin/RMout2gff3.sh sim_genome.fa.out"
echo "backbone: RepeatMasker job $jid_backbone"
fi

##########################
## 3. TEbenchtop mode 2 + mode 3, and false-positive loci
##########################
for i in {1..5};
do
res=$OUT/TEgenomeSimulator_tair10_sd_scenario1_rep${i}_result
stdgff=$res/tair10_sd_scenario1_rep${i}_repeat_annotation_out_final.gff
fai=$res/tair10_sd_scenario1_rep${i}_genome_sequence_out_final.fasta.fai
for lib in full drop25 drop50;
do
if [[ $lib == "full" ]]; then
  rmdir=$res/RepeatMasker; dep=""
else
  rmdir=$res/RepeatMasker_${lib}
  # job id from step 2, or (STEP3_ONLY) a RepeatMasker job still queued/running
  jid=${jid_rm[${i}_${lib}]:-$(squeue -h -u $user -n RMask_scenario1_rep${i}_${lib} -o %i | head -n1)}
  if [[ -n "$jid" ]]; then dep="--dependency=afterok:${jid}"; else dep=""; fi
  if [[ -z "$jid" && ! -s $rmdir/sim_genome.fa.out ]]; then echo "skip scenario1_rep${i} ${lib}: no RepeatMasker output or job"; continue; fi
fi
gname=scenario1_rep${i}
method=rm_${lib}
BM=$res/TEbenchtop_${method}
mkdir -p $BM
JOB=TEbench_${gname}_${method}

sbatch $dep <<EOF
#!/bin/bash -e
#SBATCH --job-name=${JOB}
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=$OUT/log/${JOB}_%j.out
#SBATCH --error=$OUT/log/${JOB}_%j.err

ml nextflow/24.10.1
ml apptainer/1.1
ml bedtools/2.30.0

# RepeatMasker .out -> GFF3 (writes sim_genome.fa.out.gff3 next to the .out)
cd ${rmdir}
bash ${TEBENCHTOP}/bin/RMout2gff3.sh sim_genome.fa.out
tstgff=${rmdir}/sim_genome.fa.out.gff3

# Each run gets its own launch directory (bin/ and assets/ linked from TEbenchtop), so the
# 15 concurrent nextflow runs do not share a .nextflow/ cache and lock files
LAUNCH=${BM}/launch
mkdir -p \${LAUNCH}
ln -sfn ${TEBENCHTOP}/bin \${LAUNCH}/bin
ln -sfn ${TEBENCHTOP}/assets \${LAUNCH}/assets
cd \${LAUNCH}
LOCAL_WORK=${BM}/tmp/nextflow_\${SLURM_JOB_ID}
mkdir -p \${LOCAL_WORK}

nextflow run ${TEBENCHTOP}/main.nf -profile slurm --mode 2 --fai ${fai} --gname ${gname} --method ${method} \\
    --stdgff ${stdgff} --tstgff \${tstgff} --outdir ${BM} -work-dir \${LOCAL_WORK}/m2
nextflow run ${TEBENCHTOP}/main.nf -profile slurm --mode 3 --gname ${gname} --method ${method} \\
    --stdgff ${stdgff} --tstgff \${tstgff} --outdir ${BM} -work-dir \${LOCAL_WORK}/m3

# RepeatMasker hits overlapping no simulated TE (false-positive loci)
bedtools intersect -v -a \${tstgff} -b ${stdgff} > ${BM}/fp_loci_${gname}_${method}.gff3

rm -rf \${LOCAL_WORK}
echo "${gname} ${method} done"
EOF

sleep 2
done
done
