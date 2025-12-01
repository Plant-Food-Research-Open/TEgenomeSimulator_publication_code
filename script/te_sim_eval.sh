#!/bin/bash
# Navigate to the cloned repository of TEgenomeSimulator
ml conda
conda create --name tesimeval python=3.9 
conda activate tesimeval
pip install biopython numpy scipy scikit-learn pandas matplotlib

WRK=$(pwd)
IN=$WRK/pub/input
OUT=$WRK/pub/output

############################################
# Run repeatMasker for Recovery assessment #
############################################
ml singularity/3
TOOL=../../bin/TETools
mkdir -p $OUT/TE_sim_eval
cd $OUT/TE_sim_eval
THREADS=10
RAM=15G
TIME="06:00:00"
GRP=1

# TEGSm2
prefix=tair10_ch1_m2
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
telib=$IN/athrep.updated.nonredun.noCenSatelli.fasta
JOB=RMask_$prefix
mkdir -p $OUT/TE_sim_eval/${prefix}/RepeatMasker
cd $OUT/TE_sim_eval/${prefix}/RepeatMasker
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa"

# TEGSm2+0
prefix=tair10_ch1_m0
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
telib=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
JOB=RMask_$prefix
mkdir -p $OUT/TE_sim_eval/${prefix}/RepeatMasker
cd $OUT/TE_sim_eval/${prefix}/RepeatMasker
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa"

# Garlic
prefix=Garlic
genome=$WRK/pub/output/${prefix}/tair10_garlic_sim_chr1.fa.fasta
telib=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
JOB=RMask_$prefix
mkdir -p $OUT/TE_sim_eval/${prefix}/RepeatMasker
cd $OUT/TE_sim_eval/${prefix}/RepeatMasker
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa"

################################################
# Run te_sim_eval.py & te_sim_eval_bedtools.sh #
################################################

#### TEGSm2 ####
prefix=tair10_ch1_m2
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
tegff=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_repeat_annotation_out_final.gff
rmout=$OUT/TE_sim_eval/${prefix}/RepeatMasker/sim_genome.fa.out

THREADS=1
RAM=30G
TIME="6-00:00:00"

# bedtools intersect
ml bedtools/2.30.0
out=$OUT/TE_sim_eval/${prefix}/bedtools_cut08
mkdir -p $out
cd $out
cut=0.8
$WRK/pub/script/te_sim_eval_bedtools.sh $rmout $tegff $cut

# kmer
out=$OUT/TE_sim_eval/${prefix}/kmer
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_kmer
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --outdir $out --k 6 --win 1000 --eval kmer"

# entropy
out=$OUT/TE_sim_eval/${prefix}/entropy
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_entropy
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --outdir $out --k 6 --win 1000 --eval entropy"

# compress
out=$OUT/TE_sim_eval/${prefix}/compress
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_compress
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --outdir $out --k 6 --win 1000 --eval compress"

#### TEGSm2+0 ####
prefix=tair10_ch1_m0
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
tegff=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_repeat_annotation_out_final.gff
rmout=$OUT/TE_sim_eval/${prefix}/RepeatMasker/sim_genome.fa.out

THREADS=1
RAM=30G
TIME="6-00:00:00"

for eval in kmer entropy compress
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --outdir $out --k 6 --win 1000 --eval ${eval}"
done

# bedtools intersect
ml bedtools/2.30.0
out=$OUT/TE_sim_eval/${prefix}/bedtools_cut08
mkdir -p $out
cd $out
cut=0.8
$WRK/pub/script/te_sim_eval_bedtools.sh $rmout $tegff $cut

#### Garlic ####
prefix=Garlic
genome=$WRK/pub/output/${prefix}/tair10_garlic_sim_chr1.fa.fasta
tegff=$WRK/pub/output/${prefix}/tair10_garlic_sim_chr1.fa.inserts.te.gff
rmout=$OUT/TE_sim_eval/${prefix}/RepeatMasker/sim_genome.fa.out

THREADS=1
RAM=30G
TIME="6-00:00:00"

for eval in kmer entropy compress
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --outdir $out --k 6 --win 1000 --eval ${eval}"
done

# bedtools intersect
ml bedtools/2.30.0
out=$OUT/TE_sim_eval/${prefix}/bedtools_cut08
mkdir -p $out
cd $out
sed 's/^chr1/artificial_sequence_1/g' $tegff > truth_tegff_chrmodified.gff
cut=0.8
$WRK/pub/script/te_sim_eval_bedtools.sh $rmout truth_tegff_chrmodified.gff $cut

#### Original ####
prefix=original
genome=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
tegff=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out.gff

THREADS=1
RAM=30G
TIME="6-00:00:00"

for eval in kmer entropy compress 
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --outdir $out --k 6 --win 1000 --eval ${eval}"
done


