#!/bin/bash
# Navigate to the cloned repository of TEgenomeSimulator
ml apptainer
WRK=$(pwd)
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT=$WRK/pub/output
mkdir -p $IN $OUT

# Take pseudo-chromosome only
cd $IN
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna |
grep -A1 -E '>NC_004354.4|>NT_033779.5|>NT_033778.4|>NT_037436.4|>NT_033777.3|>NC_004353.4|>NC_024512.1' > GCF_000001215.4_Release_6_chromosome_only.fna

prefix=dm6
genome=$IN/GCF_000001215.4_Release_6_chromosome_only.fna
repeat=$IN/Dmel-families.fa

JOB=TEGS_dm6
TIME="1-00:00:00"
THREADS=20
MEM=10G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"

cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
