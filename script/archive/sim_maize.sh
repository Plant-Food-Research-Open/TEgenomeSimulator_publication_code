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
WRK=$(pwd)
IN=$WRK/pub/input
OUT="../TEgenomeSimulator_pub/output"
fasta=$IN/maizeTE11122019
cd $OUT
OUT=$(pwd)

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
' "$fasta" > maizeTE11122019_uniqTEfam.fa

##############
# Rerun TEGS #
##############
ml apptainer
WRK=$(pwd)
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT="../TEgenomeSimulator_pub/output"
cd $OUT
OUT=$(pwd)
#mkdir -p $IN $OUT

prefix=maize_run3
genome=$IN/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna
repeat=$OUT/maizeTE11122019_uniqTEfam.fa

JOB=TEGS_maize_run3
TIME="3-00:00:00"
THREADS=20
MEM=50G

sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"

cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'