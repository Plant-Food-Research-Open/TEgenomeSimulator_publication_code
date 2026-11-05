cd /workspace/cflthc/script/TEgenomeSimulator_publication_code && bash script/revision/15_tair10_benchmark_repeatmasker.sh
# drop25: 88 of 353 families removed; 265 kept
# drop50: 176 of 353 families removed; 177 kept
# scenario1_rep1 drop25: RepeatMasker job 11151001
# scenario1_rep1 drop50: RepeatMasker job 11151002
# scenario1_rep2 drop25: RepeatMasker job 11151003
# scenario1_rep2 drop50: RepeatMasker job 11151004
# scenario1_rep3 drop25: RepeatMasker job 11151005
# scenario1_rep3 drop50: RepeatMasker job 11151006
# scenario1_rep4 drop25: RepeatMasker job 11151007
# scenario1_rep4 drop50: RepeatMasker job 11151008
# scenario1_rep5 drop25: RepeatMasker job 11151009
# scenario1_rep5 drop50: RepeatMasker job 11151010
# Submitted batch job 11151012
# backbone: RepeatMasker job 11151011
# Submitted batch job 11151013
# Submitted batch job 11151014
# Submitted batch job 11151015
# Submitted batch job 11151016
# Submitted batch job 11151017
# Submitted batch job 11151018
# Submitted batch job 11151019
# Submitted batch job 11151020
# Submitted batch job 11151022
# Submitted batch job 11151023
# Submitted batch job 11151024
# Submitted batch job 11151025
# Submitted batch job 11151026
# Submitted batch job 11151027
# Submitted batch job 11151028

# Resubmit the scoring step:
cd /workspace/cflthc/script/TEgenomeSimulator_publication_code && STEP3_ONLY=1 bash script/revision/15_tair10_benchmark_repeatmasker.sh
# Submitted batch job 11157034
# Submitted batch job 11157037
# Submitted batch job 11157041
# Submitted batch job 11157045
# Submitted batch job 11157048
# Submitted batch job 11157052
# Submitted batch job 11157056
# Submitted batch job 11157060
# Submitted batch job 11157063
# Submitted batch job 11157070
# Submitted batch job 11157074
# Submitted batch job 11157081
# Submitted batch job 11157088
# Submitted batch job 11157092
# Submitted batch job 11157098

# Knit Rmd
sbatch --cpus-per-task=1 --mem=32G -t 04:00:00 -J bm_knit -o bm_knit.%j.out -e bm_knit.%j.err --wrap "cd /workspace/cflthc/script/TEgenomeSimulator_publication_code/script/revision && Rscript -e \"rmarkdown::render('25.1_item23_28_benchmark_repeatmasker.Rmd')\""
Submitted batch job 11158334