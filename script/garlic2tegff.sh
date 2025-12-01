#!/bin/bash

awk -F'\t' '
function clamp(v, lo, hi){ return v<lo?lo:(v>hi?hi:v) }

NF > 0 {
  start0 = $1
  end    = $2
  idx    = $3
  info   = $4

  sub(/\[.*/,"", info)

  n = split(info, parts, /;/)
  left = parts[1]
  right = (n>1 ? parts[2] : "")

  split(left, L, /:/)
  TEfam   = L[1]
  TEclass = L[2]
  orient  = L[3]
  pctdiv  = (L[4]  ? L[4]  : 0) + 0
  pctins  = (L[5]  ? L[5]  : 0) + 0
  pctdel  = (L[6]  ? L[6]  : 0) + 0
  fraglen = (L[7]  ? L[7]  : 0) + 0
  brk     = (L[8]  ? L[8]  : 0) + 0

  split(right, R, /:/)
  inslvl  = (R[1] ? R[1] : 0) + 0
  p_tran  = (R[2] ? R[2] : 0) + 0
  p_trv   = (R[3] ? R[3] : 0) + 0
  p_del   = (R[4] ? R[4] : 0) + 0
  p_ins   = (R[5] ? R[5] : 0) + 0

  broad = TEclass
  sub(/\/.*/,"", broad)

  identity_pct = 100 - (p_tran + p_trv + p_trv + p_ins)
  identity_pct = clamp(identity_pct, 0, 100)
  identity     = identity_pct / 100.0

  denom = 1.0 - (p_del/100.0) + (p_ins/100.0)
  integrity = (denom > 0 ? ( (end - start0 + 1) / denom ) : 0)
  integrity_len_norm = (end - start0 + 1) > 0 ? integrity / (end - start0 + 1) : 0
  integrity_len_norm = clamp(integrity_len_norm, 0, 1)

  chr     = "chr1"
  source  = "Garlic"
  feature = broad
  start1  = start0 + 1
  score   = "."
  strand  = orient
  phase   = "."

  id      = "TE" idx "_" TEfam
  attrs   = "ID=" id ";Classification=" TEclass ";Identity=" sprintf("%.4f", identity) ";Integrity=" sprintf("%.4f", integrity_len_norm)

  OFS = "\t"
  print chr, source, feature, start1, end, score, strand, phase, attrs
}
' "$1"
