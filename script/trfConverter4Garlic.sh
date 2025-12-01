#!/bin/bash
# usage: trfConverter4Garlic.sh trf_out.dat

input=$1

awk '
  /^Sequence:/ {
    # Extract the sequence name (first word after "Sequence:")
    split($2, seq_id, " ")
    current_seq = seq_id[1]
  }
  /^[0-9]/ {
    # Data line: starts with a number
    start = $1
    end = $2
    period = $3
    matched = $6
    indel = $7
    consensus = $NF
    printf "%s %s %s trf . %s . %s %s %s\n", current_seq, start, end, period, matched, indel, consensus
  }
' "$input"
