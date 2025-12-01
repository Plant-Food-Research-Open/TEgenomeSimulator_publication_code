#!/usr/bin/env bash
# Use RepeatMasker.out and repeat's FAI file to generate informative gff3 file. 
DATE=$(date -u)

in_rmout=$1
in_fai=$2
out_gff=$3
col2label=$4

# Extract start, end, TE family, classification, and identity from RepeatMasker.out
awk -v c2=$col2label 'NR>3, idt=(100-$2)/100 {print $5"\t"c2"\t"$11"\t"$6"\t"$7"\t"$1"\t"$9"\t""\.""\t""ID="$10"_"$5"_"$6"-"$7";Name="$10";Classification="$11";Identity="idt}' $in_rmout > temp1.gff3

# Add hashtag info
awk 'BEGIN{FS=OFS="\t"} {gsub("C", "-", $7)} 1' temp1.gff3 |\
sed '1i ##Column 6: Smith-Waterman score extracted from the first column of RepeatMasker.out' |\
sed '1i ##Integrity: Sequence integrity (0-1) between the library sequence and the target region.' |\
sed '1i ##Identity: Sequence identity (0-1) between the library sequence and the target region.' |\
sed "1i ##date: $DATE" |\
sed '1i ##gff-version 3' > temp2.gff3

# -------------------------------
# Process GFF file
# -------------------------------
# The curly braces { ... } group multiple commands together.
# At the end we write: } > "$out_file"
# That means the combined output of ALL commands inside { ... }
# will be redirected once into $out_file.
# Without the braces, only the last command would redirect, overwriting the file.
{
    # Step 1: print header lines (those starting with '#') unchanged
    grep "^#" "temp2.gff3"

    # Step 2: process non-header lines through awk
    grep -v "^#" "temp2.gff3" | awk -v OFS="\t" -v fai_file="$in_fai" '
        BEGIN {
            # Load TE family lengths from the fai file into an array
            while ((getline < fai_file) > 0) {
                split($1, a, "#");         # separate name at "#"
                fam = a[1];                # TE family = part before "#"
                fai_len[fam] = $2;         # store length
            }
            close(fai_file);
        }
        {
            start = $4; 
            end   = $5;
            te_length = end - start + 1;   # TE instance length

            # Extract "Name=" attribute from column 9
            if (match($9, /Name=([^;]+)/, arr)) {
                te_family = arr[1];        # capture TE family name
                fam_len = fai_len[te_family];

                if (fam_len > 0) {
                    # Calculate integrity (rounded to 3 decimals)
                    integrity = sprintf("%.3f", te_length / fam_len);
                    $9 = $9 ";Integrity=" integrity;
                } else {
                    # If family not found in fai, mark as NA
                    $9 = $9 ";Integrity=NA";
                }
            }
            print $0;  # print the modified GFF line
        }
    '
# Close the command group. Apply output redirection to the entire block.
} > "$out_gff"

rm temp1.gff3 temp2.gff3