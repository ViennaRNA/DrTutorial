#!/bin/bash
#

SEQUENCE_DIR="sequences"
BACKBONE_COL="0.7 0.6 0.6"
PAIR_COL="0.9 0.6 0.4"
NATIVE_STEM_COL="0.4 0.6 0.4"
NATIVE_COL="0.7 0.9 0.7"
ANNOT_COL="0. 0.3 0.6"
POS_COL="0.5 0.5 0.5"
LABEL_FONT="0.1 0.3 0.1 setrgbcolor (Times) 2 LabelFont"

ANNOT_H1_STEM="3 25 ${NATIVE_COL} Fomark"
ANNOT_H1_CONST="6 22 8 20 ${NATIVE_STEM_COL} BFmark"
ANNOT_H1A_STEM="4 22 ${NATIVE_COL} Fomark"
ANNOT_H1A_CONST="8 18 10 16 ${NATIVE_STEM_COL} BFmark"
ANNOT_H1B_STEM="4 23 ${NATIVE_COL} Fomark"
ANNOT_H1B_CONST="6 21 8 19 ${NATIVE_STEM_COL} BFmark" # 8 31 9 30 ${NATIVE_STEM_COL} BFmark 12 27 14 25 ${NATIVE_STEM_COL} BFmark 16 23 17 22 ${NATIVE_STEM_COL} BFmark"
ANNOT_H1C_STEM="3 38 ${NATIVE_COL} Fomark"
ANNOT_H1C_CONST="8 31 9 30 ${NATIVE_STEM_COL} BFmark 12 27 14 25 ${NATIVE_STEM_COL} BFmark 16 23 17 22 ${NATIVE_STEM_COL} BFmark"
ANNOT_H2_STEM="46 69 ${NATIVE_COL} Fomark"
ANNOT_H2_CONST="46 69 48 67 ${NATIVE_STEM_COL} BFmark 53 62 55 60 ${NATIVE_STEM_COL} BFmark"
ANNOT_S1_CONST="38 74 41 71 ${NATIVE_STEM_COL} BFmark"
ANNOT_S2_CONST="28 85 31 82 ${NATIVE_STEM_COL} BFmark"
ANNOT_S3_CONST="22 93 24 91 ${NATIVE_STEM_COL} BFmark"
ANNOT_S4_CONST="4 113 19 98 ${NATIVE_STEM_COL} BFmark"
ANNOT_H2A_STEM="36 57 ${NATIVE_COL} Fomark 36 57 38 55 ${NATIVE_STEM_COL} BFmark 42 53 44 51 ${NATIVE_STEM_COL} BFmark"
ANNOT_H3A_STEM="72 90 ${NATIVE_COL} Fomark 72 90 75 87 ${NATIVE_STEM_COL} BFmark 77 85 79 83 ${NATIVE_STEM_COL} BFmark"
ANNOT_H3_STEM="87 105 ${NATIVE_COL} Fomark"
ANNOT_H3_CONST="87 105 90 102 ${NATIVE_STEM_COL} BFmark"

ANNOT_SETTINGS="${ANNOT_COL} setrgbcolor 2 setlinewidth"

RNAplot \
  -t4 \
  --pre="${ANNOT_H1_STEM} ${ANNOT_H1_CONST} gsave ${ANNOT_SETTINGS} 21 cmark 21 1 -0.4 (U21C) Label 22 cmark 22 1 -0.4 (C22U) Label ${POS_COL} setrgbcolor (Helvetica) 0.8 LabelFont 10 -1.5 -0.4 (10) Label 20 0.5 -0.4 (20) Label 1 -1.5 -0.4 (5') Label 26 1 -0.4 (3') Label ${LABEL_FONT} 6 -3.5 -0.5 (H1) Label grestore" \
  ${SEQUENCE_DIR}/SRPn-H1.fa

RNAplot \
  -t4 \
  --pre="${ANNOT_H1A_STEM} ${ANNOT_H1A_CONST} gsave ${ANNOT_SETTINGS} 21 cmark 21 1 -0.4 (U21C) Label 22 cmark 22 -1 -1.5 (C22U) Label ${POS_COL} setrgbcolor (Helvetica) 0.8 LabelFont 10 -1.5 -0.4 (10) Label 20 0.5 -0.4 (20) Label 1 -1.5 -0.4 (5') Label 24 1 -0.4 (3') Label ${LABEL_FONT} 6 -4.5 -0.5 (H1a) Label grestore" \
  ${SEQUENCE_DIR}/SRPn-H1a.fa

RNAplot \
  -t4 \
  --pre="${ANNOT_H1B_STEM} ${ANNOT_H1B_CONST} gsave ${ANNOT_SETTINGS} 21 cmark 21 1 -0.4 (U21C) Label 22 cmark 22 1 -0.4 (C22U) Label ${POS_COL} setrgbcolor (Helvetica) 0.8 LabelFont 10 -1.5 -0.4 (10) Label 20 0.5 -0.4 (20) Label 1 -1.5 -0.4 (5') Label 25 1 -0.4 (3') Label ${LABEL_FONT} 6 -4.5 -0.5 (H1b) Label grestore" \
  ${SEQUENCE_DIR}/SRPn-H1b.fa

RNAplot \
  -t4 \
  --pre="${ANNOT_H1C_STEM} ${ANNOT_H1C_CONST} gsave ${ANNOT_SETTINGS} 35 cmark 37 cmark 35 1 -0.4 (U35C) Label 37 1 -0.4 (U37C) Label ${POS_COL} setrgbcolor (Helvetica) 0.8 LabelFont 10 -0.5 -1 (10) Label 20 -0.5 0.5 (20) Label 30 0.5 0 (30) Label 1 -1.5 -0.4 (5') Label 40 1 -0.4 (3') Label ${LABEL_FONT} 11 -4.5 -0.5 (H1c) Label  grestore" \
  ${SEQUENCE_DIR}/SRPn-H1c.fa

RNAplot \
  -t4 \
  --pre="${ANNOT_H1_STEM} ${ANNOT_H1_CONST} ${ANNOT_H2_STEM} ${ANNOT_H2_CONST} ${ANNOT_S1_CONST} ${ANNOT_S2_CONST} ${ANNOT_H3_STEM} ${ANNOT_H3_CONST} gsave ${ANNOT_COL} setrgbcolor 2 setlinewidth 21 cmark 22 cmark 35 cmark 37 cmark 93 cmark ${POS_COL} setrgbcolor (Helvetica) 0.8 LabelFont 10 -1.5 -0.4 (10) Label 20 0.5 -0.4 (20) Label 30 -1.5 -0.4 (30) Label 40 -1.5 -0.4 (40) Label 50 -1.5 -0.4 (50) Label 60 0.5 -0.5 (60) Label 70 0.5 -0.4 (70) Label 80 0.5 -0.4 (80) Label 90 -1.5 -0.4 (90) Label 100 0.5 -0.4 (100) Label 1 -1.5 -0.4 (5') Label 106 1 -0.4 (3') Label ${LABEL_FONT} 6 -3.5 -0.5 (H1) Label 73 1 -0.5 (S1) Label 83 1 -0.5 (S2) Label 50 -4.5 -0.5 (H2) Label 103 1 -0.5 (H3) Label grestore" \
  ${SEQUENCE_DIR}/SRPn-H1-H2-S1-S2-H3.fa

RNAplot \
  -t4 \
  --pre="${ANNOT_H1B_STEM} ${ANNOT_H1B_CONST} ${ANNOT_H2A_STEM} ${ANNOT_H3A_STEM} gsave ${ANNOT_COL} setrgbcolor 2 setlinewidth 21 cmark 22 cmark 35 cmark 37 cmark 93 cmark ${POS_COL} setrgbcolor (Helvetica) 0.8 LabelFont 10 -1.5 -0.4 (10) Label 20 0.5 -0.4 (20) Label 30 -0.4 -1.2 (30) Label 40 -1.3 -0.4 (40) Label 50 0.5 -0.4 (50) Label 60 0.5 -0.1 (60) Label 70 -1.3 -0.4 (70) Label 80 -0.6 0.5 (80) Label 90 0.5 -0.4 (90) Label 1 -1.5 -0.4 (5') Label 95 1 -0.4 (3') Label ${LABEL_FONT} 6 -4.5 -0.5 (H1b) Label 55 1 -0.5 (H2a) Label 88 1 -0.5 (H3a) Label grestore" \
  ${SEQUENCE_DIR}/SRPn-H1b-H2a-H3a.fa

RNAplot \
  -t4 \
  --pre="${ANNOT_H2_STEM} ${ANNOT_H2_CONST} ${ANNOT_S1_CONST} ${ANNOT_S2_CONST} ${ANNOT_S3_CONST} ${ANNOT_S4_CONST} gsave ${ANNOT_COL} setrgbcolor 2 setlinewidth 21 cmark 21 -3.5 -0.4 (U21C) Label 22 cmark 22 -3.5 -0.4 (C22U) Label 35 cmark 37 cmark 35 -3.5 -0.4 (U35C) Label 37 -3.5 -0.4 (U37C) Label 93 cmark 93 1 -0.4 (G93A) Label ${POS_COL} setrgbcolor (Helvetica) 0.8 LabelFont 10 -1.5 -0.4 (10) Label 20 -1.5 -0.4 (20) Label 30 -1.5 -0.4 (30) Label 40 -1.5 -0.4 (40) Label 50 -1.5 -0.4 (50) Label 60 0 -1.2 (60) Label 70 0.5 -0.4 (70) Label 80 0.5 -0.4 (80) Label 90 0.5 -0.4 (90) Label 100 0.5 -0.4 (100) Label 110 0.5 -0.4 (110) Label 1 -1.5 -0.4 (5') Label 117 1 -0.4 (3') Label ${LABEL_FONT} 50 -4.5 -0.5 (H2) Label 73 2 -0.5 (S1) Label 84 2 -0.5 (S2) Label 92 2 -0.5 (S3) Label 106 2 -0.5 (S4) Label grestore" \
  ${SEQUENCE_DIR}/SRPn-H2-S1-S2-S3-S4.fa

# rename files
for f in \
    SRPn-H1_ss.ps \
    SRPn-H1a_ss.ps \
    SRPn-H1b_ss.ps \
    SRPn-H1c_ss.ps \
    SRPn-H1-H2-S1-S2-H3_ss.ps \
    SRPn-H1b-H2a-H3a_ss.ps \
    SRPn-H2-S1-S2-S3-S4_ss.ps
do
  # post-process
  sed -i 's#/fsize  14#/fsize  18#g' $f
  sed -i 's#/outlinecolor {0.2 setgray}#/outlinecolor { '"${BACKBONE_COL}"' setrgbcolor}#g' $f
  sed -i 's#/paircolor    {0.2 setgray} bind def#/paircolor     {'"${PAIR_COL}"' setrgbcolor} bind def#g' $f
  sed -i 's#  0.7 setlinewidth#  5 setlinewidth#g' $f
  sed -i 's#/seqcolor     {0   setgray}#/seqcolor     {0.3   setgray}#' $f

  # increase plot area
  sed -i 's#%%BoundingBox: 0 0 700 700#%%BoundingBox: -70 -70 770 770#' $f

#  sed -i 's#^drawoutline#%drawoutline#g' $f

  /usr/bin/mv $f `basename $f _ss.ps`.eps
done
