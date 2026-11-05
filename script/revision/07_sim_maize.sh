#!/bin/bash
# Navigate to the cloned repository of TEgenomeSimulator

########################################
# Found out that the curated maize TE family acquired from EDTA GitHub have duplicated TE family name.
# E.g:
# >ZM00010_consensus#DNA/hAT
# >ZM00010_consensus#DNA/CACTA
# >ZM00010_consensus#DNA/Harbinger
# >ZM00010_consensus#DNA/Tc1-Mariner
# This cause problem for simulation
# Need to modify the TE fasta file to get around with it.
########################################
# Navigate to the cloned repository of "TEgenomeSimulator_publication_code"
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # remove this line before git push
cd $WRK
IN=$WRK/input
OUT=$WRK/output
fasta=$IN/maizeTE11122019

# 1. Count total headers
total=$(grep -c '^>' "$fasta")

# 2. Determine how many digits are needed
digits=$(awk -v t="$total" 'BEGIN{
  d=length(t);
  if (d<2) d=2;  # minimum 2 digits
  print d
}')

# 3. Add numbers before first '#'
awk -v digits="$digits" '
BEGIN {i=0}
# Header lines
/^>/ {
    i++
    header=$0
    # Find position of first #
    pos=index(header, "#")
    if (pos==0) {
        # no # found, just append at the end of the line
        printf "%s_%0*d\n", header, digits, i
    } else {
        # split at #
        before=substr(header,1,pos-1)
        after=substr(header,pos)
        printf "%s_%0*d%s\n", before, digits, i, after
    }
    next
}
# Sequence lines unchanged
{print}
' "$fasta" > $IN/maizeTE11122019_uniqTEfam.fa

##############
# Rerun TEGS #
##############
# GITHUB_REPO=cloned repository of "TEgenomeSimulator"
# Navigate to the cloned repository of "TEgenomeSimulator_publication_code"
ml apptainer
ml samtools/1.20
#TEgenomeSimulator=/workspace/${user}/script/TEgenomeSimulator/TEgenomeSimulator_dev.sif # remove this line before git push
WRK=/workspace/${user}/script/TEgenomeSimulator_publication_code # remove this line before git push
cd $WRK
#TEgenomeSimulator=$GITHUB_REPO/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=maize
genome=$IN/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna
repeat=$IN/maizeTE11122019_uniqTEfam.fa

JOB=TEGS_${prefix}
TIME="3-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"  # remove this line before git push
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $OUT/log/$JOB.out -e $OUT/log/$JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
# Submitted batch job 10737365

cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'