#!/bin/bash
#!/bin/bash

# Inputs
rm_out=$1
truth_gff=$2
cutoff=$3 # minimum reciprocal overlap

# 1. Convert RepeatMasker .out to BED, ensure start < end
awk 'NR>3 {
    start=$6-1; end=$7;
    if (start > end) {tmp=start; start=end; end=tmp}
    print $5 "\t" start "\t" end "\t" $10
}' "$rm_out" > repeatmasker.bed

# 2. Convert truth GFF to BED, ensure start < end
awk '{
    start=$4-1; end=$5;
    if (start > end) {tmp=start; start=end; end=tmp}
    print $1 "\t" start "\t" end "\t" $9
}' "$truth_gff" > truth.bed

# 3. Intersect: with overlap and without overlap
bedtools intersect -u -f $cutoff -r -a truth.bed -b repeatmasker.bed > truth_with_overlap.bed
bedtools intersect -v -f $cutoff -r -a truth.bed -b repeatmasker.bed > truth_without_overlap.bed


# 4. Extract integrity values
extract_integrity () {
  awk -v cat=$1 '{
    integ=""; ident="";
    if ($4 ~ /Integrity=/) {
      match($4,/Integrity=([0-9.]+)/,a); integ=a[1];
    }
    if ($4 ~ /Identity=/) {
      match($4,/Identity=([0-9.]+)/,b); ident=b[1];
    }
    if (integ != "" && ident != "") {
      print cat, integ, ident;
    }
  }' $2
}

    # Run for both categories
echo -e "category\tintegrity\tidentity" > with_overlap_integrity.tsv
extract_integrity with_overlap truth_with_overlap.bed >> with_overlap_integrity.tsv

echo -e "category\tintegrity\tidentity" > without_overlap_integrity.tsv
extract_integrity without_overlap truth_without_overlap.bed >> without_overlap_integrity.tsv

# 5. Summarize per bin with labeled intervals
cat with_overlap_integrity.tsv without_overlap_integrity.tsv | \
sed -e 's/^category.*//g' -e '/^$/d' | \
  awk '
  NR==1 && $1=="category" {next}                                                   # Skip header lines
  {
    integ_bin = int($2*10);                                                        # Convert integrity (col 2) into a bin index 0-9
    if (integ_bin==10) integ_bin=9;                                                # Ensure values == 1.0 go into bin 9
    div_bin   = int($3*10);                                                        # Convert identity (col 3) into a bin index 0-9
    if (div_bin==10) div_bin=9;                                                    # Same for identity = 1.0

    key = integ_bin "_" div_bin;                                                   # Create 2D key combining both bin indices
    total[key]++;                                                                  # Count all records in this bin
    if ($1=="with_overlap") recovered[key]++;                                      # Count only "with_overlap" ones
  }
  END{
    print "integrity_interval","identity_interval","total","recovered","rate";     
    for (i=0;i<10;i++) {                                                           # Loop over integrity bins
      low_i = i/10; high_i = (i+1)/10;                                             # Compute bin boundaries
      if (i<9) integ_label = sprintf("[%.1f-%.1f)", low_i, high_i);                # Inclusive–exclusive
      else     integ_label = sprintf("[%.1f-%.1f]", low_i, high_i);                # Last bin inclusive
      for (j=0;j<10;j++) {                                                         # Loop over identity bins
        low_j = j/10; high_j = (j+1)/10;
        if (j<9) div_label = sprintf("[%.1f-%.1f)", low_j, high_j);
        else     div_label = sprintf("[%.1f-%.1f]", low_j, high_j);
        key = i "_" j;                                                             # Rebuild key
        tot = (key in total ? total[key] : 0);                                     # Total entries in this bin
        rec = (key in recovered ? recovered[key] : 0);                             # Recovered entries
        rate = (tot>0 ? rec/tot : 0);                                              # Recovery rate (avoid division by zero)
        print integ_label, div_label, tot, rec, rate;
      }
    }
  }' OFS="\t" > recovery_by_bin2D.tsv

# 6. Calculate numbers
total=$(wc -l < truth.bed)
with_overlap=$(wc -l < truth_with_overlap.bed)
no_overlap=$(( total - with_overlap ))
ratio=$(echo "$no_overlap / $total" | bc -l)

# 7. Print tab-delimited results
echo -e "Total_truth\tTruth_with_overlap\tTruth_without_overlap\tRatio_non_overlap" > overlap_eval_stat.txt
echo -e "$total\t$with_overlap\t$no_overlap\t$ratio" >> overlap_eval_stat.txt

