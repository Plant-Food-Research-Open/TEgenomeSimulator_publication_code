#!/bin/bash

# genome
# Extract chr1 of TAIR10
WRK=$(pwd)
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna |
grep -A1 -E 'chromosome 1' > $WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna

# extract chr1 TE annotation for visulisation
cd $WRK/pub/input
fullRMout=../output/TEgenomeSimulator_tair10_ch1_m2_result/repeatmasker/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out
cp $fullRMout .
    # convert .out to .gff
in_rmout=GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out
in_fai=athrep.updated.nonredun.noCenSatelli.fasta.stitched.fai
out_gff=GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out.gff
col2label="RepeatMasker"
$WRK/pub/script/RMout2gff3v2.sh $in_rmout $in_fai $out_gff $col2label

# repeat
cp ../../../../scratch/TEGE/TEgenomeSimulator/tair10/output/TEgenomeSimulator_tair10_run6_result/repeatmasker/GCF_000001735.4_TAIR10.1_genomic.fna.align ./

# TRF.out
ml trf/4.07
cd $WRK/pub/input
trf GCF_000001735.4_TAIR10.1_genomic_chr1.fna 2 5 7 80 10 50 2000 -h -d | trf -ngs -

JOB=TRF_tair10
TIME="3-00:00:00"
THREADS=1
MEM=15G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "trf $WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna 2 5 7 80 10 50 2000 -h -d | trf -ngs -"

# Convert TRF output .dat file to the format required for Garlic
../script/trfConverter4Garlic.sh GCF_000001735.4_TAIR10.1_genomic.fna.2.5.7.80.10.50.2000.dat > GCF_000001735.4_TAIR10.1_genomic.fna.2.5.7.80.10.50.2000.dat.reformated

# Gene table in UCSC ensGenes.txt format
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.61.gff3.gz
gzip -d -k Arabidopsis_thaliana.TAIR10.61.gff3.gz
JOB=gff32ensGene_tair10_chr1
TIME="2-00:00:00"
THREADS=1
MEM=10G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python ../../../../bin/gff32ensGene.py -i Arabidopsis_thaliana.TAIR10.61.gff3"

# Check the length of tair10 chr1
ml samtools
cd $WRK/pub/input/
samtools faidx $WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna 

# Convert the TE lib I used to embl format and replace the RepBase db
infa=../../../../scratch/TEGE/TEgenomeSimulator/tair10/output/TEgenomeSimulator_tair10_run6_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
cp $infa ./
outembl=$WRK/pub/input/athrep.updated.nonredun.noCenSatelli.fasta.stitched.embl
python3 - <<EOF
from Bio import SeqIO
records = []
for r in SeqIO.parse("$infa", "fasta"):
    r.annotations["molecule_type"] = "DNA"
    records.append(r)
SeqIO.write(records, "$outembl", "embl")
EOF

# The converted TE lib in embl format has to be stored under the folder "RepBase" and named as "RepeatMaskerLib.embl"
OUT=$WRK/pub/output/Garlic
cd $OUT
rm -rf RepBase/*
cp $outembl RepBase/
mv RepBase/athrep.updated.nonredun.noCenSatelli.fasta.stitched.embl RepBase/RepeatMaskerLib.embl
# Replacing '#' to ' ' otherwise no TE insertion recorded by Garlic
sed 's/#/ /g' ./RepBase/RepeatMaskerLib.embl > ./RepBase/temp.embl
mv ./RepBase/temp.embl ./RepBase/RepeatMaskerLib.embl

# Run createModel.pl
# Use RepeatMasker's '.align' file instead of '.out' file
OUT=$WRK/pub/output/Garlic
mkdir -p $OUT
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
genome=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
RMout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna.align
TRFout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna.2.5.7.80.10.50.2000.dat.reformated
ensGene=$WRK/pub/input/ensGene_table.tsv
perl $GarlicDir/createModel.pl -m tair10 -d $OUT -f $genome -r $RMout -t $TRFout -g $ensGene

# Run createFakeSequence.pl
# The embl file converted from A. thaliana TE lib was used as database, instead of the one downloaded from repbase.
JOB=GarlicSim_tair10
TIME="6-00:00:00"
THREADS=1
MEM=20G
OUT=$WRK/pub/output/Garlic
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
size=$(awk '{print $2}' ../../input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.fai) # echo $size # 30427671
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "perl $GarlicDir/createFakeSequence.pl -m tair10 -s $size -d . -o tair10_garlic_sim_chr1.fa -v 2> createFakeSequence.log"


ml samtools
samtools faidx tair10_garlic_sim_chr1.fa.fasta

# The records of TE and Simple Repeat insertions in the .inserts file are in different format. Need to saperate them before importing them to R.
# Note: Garlic analyses CDS profile but doesn't simulate CDS in the synthetic genome.
WRK=$(pwd)
OUT=$WRK/pub/output/Garlic
cd $OUT
grep -v "#" tair10_garlic_sim_chr1.fa.inserts | grep -v "SIMPLE" > tair10_garlic_sim_chr1.fa.inserts.te
grep "SIMPLE" tair10_garlic_sim_chr1.fa.inserts > tair10_garlic_sim_chr1.fa.inserts.simple
grep -v "#" tair10_2_garlic_sim_chr1.fa.inserts | grep -v "SIMPLE" > tair10_2_garlic_sim_chr1.fa.inserts.te
grep "SIMPLE" tair10_2_garlic_sim_chr1.fa.inserts > tair10_2_garlic_sim_chr1.fa.inserts.simple

# Convert .insert output to gff format
cd $OUT
../../script/garlic2tegff.sh tair10_garlic_sim_chr1.fa.inserts.te > tair10_garlic_sim_chr1.fa.inserts.te.gff
../../script/garlic2simprepgff.sh tair10_garlic_sim_chr1.fa.inserts.simple > tair10_garlic_sim_chr1.fa.inserts.simple.gff

# The .inserts output looks like:
### ARTIFICIAL SEQUENCE 1 ###
#(INI   END]-zero_based NUM     REPEAT  REPEAT_EVOL
#20462816        20463231        1       VANDAL3:DNA/MuDR:+:21.88:4.65:3.88:411:1;1:9.00:10.22:4.62:3.65[$inserted_TE_sequence] # Note the coordinates are zero-based, so length = END - INI. 
# The followings are information from the script `createFakeSequence.pl` (from line 1156)
#push @inserts,
      #    "$pos\t$pos_end\t$urep\t$new\[$seq\]\t$info\n$con\n$mat\n$mut\n";
      # RMH: The transitions/transversion/insertions/deletions are now calcualted
      #      by evolveRepeat and reported after the canonical repeat pattern.
      #      The format is:
      #         prototype_details;instance_details[sequence],...
      #      Where prototype details are:
      #         family:class:orientation:%div:%ins:%del:fraglen:break
      #
      #         family: The family identifier for the TE
      #         class: The TE classification
      #         orientation: '+','-'
      #         %div: The percent divergence of the prototype being simulated
      #         %ins: The percent insertion of the prorotype being simulated
      #         %del: The percent deletion of the prototype being simulated
      #         fraglen: The length of the TE fragment being simulated           # Note: this is the length before adding small insertion and deletion!
      #         break: The number of fragments for this prototype
      #
      #      and instance details are:
      #         level:%transitions:%transversions:%insertions:%deletions
      #
      #         level: For fragmented insertions this indicates the level of 
      #                the fragment (starting from 1).  e.g An insertion of a
      #                MLT1 in an AluSx would have three fragments listed as:
      #                AluSx:::::::;1::::
      #                MLT11:::::::;2::::
      #                AluSx:::::::;1::::
      #
      #         %transitions: The percent transitions (relative to the consensus
      #                       fragment length).
      #         %transversions: The percent transversions (relative to the consensus
      #                       fragment length).
      #         %deletions: The percent deletions (relative to the consensus
      #                       fragment length).
      #         %insertions: The percent insertions (relative to the consensus
      #                       fragment length).
      # e.g:
      # AluSz:SINE/Alu:-:16.4:0.0:0.6:308:1;1:7.14:8.12:0.00:0.32[GGCCGGGGGCGG...]
      #    
# 