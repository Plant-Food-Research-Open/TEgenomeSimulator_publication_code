#!/bin/bash

##########################
## mean identity 90-100 ##
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_90_100
minidn=90
maxidn=100
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##########################
## mean identity 80-90 ##
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_80_90
minidn=80
maxidn=90
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=1
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##########################
## mean identity 70-80 ##
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_70_80
minidn=70
maxidn=80
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=1
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##########################
## mean identity 85-100 ##
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_85_100
minidn=85
maxidn=100
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 80-95 ##
#########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_80_95
minidn=80
maxidn=95
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 75-90 ##
#########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_75_90
minidn=75
maxidn=90
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 70-85 ##
#########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_70_85
minidn=70
maxidn=85
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## mean identity 80-100 ##
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_80_100
minidn=80
maxidn=100
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 75-95 ##
#########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_75_95
minidn=75
maxidn=95
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 70-90 ##
#########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_70_90
minidn=70
maxidn=90
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=1
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 70-95 ##
#########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_70_95
minidn=70
maxidn=95
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=1
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##########################
## mean identity 70-100 ##
##########################
ml singularity/3
ml samtools/1.20
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
WRK=$(pwd)
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_70_100
minidn=70
maxidn=100
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=1
MEM=50G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT"

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'