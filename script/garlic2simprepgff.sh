#!/usr/bin/env bash
# Usage: ./garlic2simprepgff.sh input.txt > output.gff

awk -F'\t' '
NF >= 4 {
    start0 = $1 + 0
    end    = $2 + 0
    idx    = $3
    info   = $4

    # Split the info field: SIMPLE:pattern:orientation:length:...
    split(info, arr, ":")
    type       = arr[1]       # SIMPLE
    orientation = arr[3]      # + or -

    chr    = "chr1"
    source = "Garlic"
    feature= type
    start1 = start0 + 1
    score  = "."
    strand = (orientation == "" ? "." : orientation)
    phase  = "."

    # Build attributes
    attrs = "ID=" type idx ";Classification=" type

    OFS = "\t"
    print chr, source, feature, start1, end, score, strand, phase, attrs
}
' "$1"
