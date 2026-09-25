#!/bin/bash
#BLAT_primers.sh

echo "Hi, $USER. Compare BCs to original dataset."
echo

#BLAT required before proceeding

#use Conda environment seqanalysis

#convert all necessary files to 2bit format

DATADIR=seq_files
AMPPRIMERDIR=ampprimers-skpp15
FILES="$DATADIR/*.fasta"

for f in $FILES
do
  #echo "Processing $f file..."
  result_string="${f/fasta/2bit}"
  faToTwoBit $f $result_string
done

#load full 25mer set
#Align the skpp20 FWD 20-mers.
blat seq_files/bc25mer.240k.2bit seq_files/skpp20-forward.fasta BCs_psl/skpp20FWD.psl -tileSize=18 -stepSize=1 -minMatch=1 -minScore=0 -minIdentity=100 -noHead -t=dna -q=dna

sleep 8

pslPretty BCs_psl/skpp20FWD.psl seq_files/bc25mer.240k.2bit seq_files/skpp20-forward.fasta BCs_psl/skpp20FWD_pretty.txt

#Align the skpp20 REV 20-mers.
blat seq_files/bc25mer.240k.2bit seq_files/skpp20-reverse.fasta BCs_psl/skpp20REV.psl -tileSize=18 -stepSize=1 -minMatch=1 -minScore=0 -minIdentity=100 -noHead -t=dna -q=dna

sleep 8

pslPretty BCs_psl/skpp20REV.psl seq_files/bc25mer.240k.2bit seq_files/skpp20-reverse.fasta BCs_psl/skpp20REV_pretty.txt

##############################################

#Align the skpp15 15-mers FWD
blat seq_files/bc25mer.240k.2bit seq_files/skpp15-forward.fasta BCs_psl/skpp15FWD.psl -tileSize=15 -stepSize=1 -minMatch=1 -minScore=0 -minIdentity=100 -noHead -t=dna -q=dna

sleep 8

pslPretty BCs_psl/skpp15FWD.psl seq_files/bc25mer.240k.2bit seq_files/skpp15-forward.fasta BCs_psl/skpp15FWD_pretty.txt

#Align the skpp15 15-mers REV
blat seq_files/bc25mer.240k.2bit seq_files/skpp15-reverse.fasta BCs_psl/skpp15REV.psl -tileSize=15 -stepSize=1 -minMatch=1 -minScore=0 -minIdentity=100 -noHead -t=dna -q=dna

sleep 8

pslPretty BCs_psl/skpp15REV.psl seq_files/bc25mer.240k.2bit seq_files/skpp15-reverse.fasta BCs_psl/skpp15REV_pretty.txt


##############################################

#Align 384 BC set.
blat seq_files/bc25mer.240k.2bit seq_files/filt_prim_12nt_Lev_3_Tm_40_42_GC_45_55_SD_2_mod_restriction_trim.fasta BCs_psl/bc384.psl -tileSize=12 -stepSize=1 -minMatch=1 -minScore=0 -minIdentity=100 -noHead -t=dna -q=dna

sleep 8

pslPretty BCs_psl/bc384.psl seq_files/bc25mer.240k.2bit seq_files/filt_prim_12nt_Lev_3_Tm_40_42_GC_45_55_SD_2_mod_restriction_trim.fasta BCs_psl/bc384_pretty.txt

##############################################

#Align 1536 BC set.
#bbmap.sh t=4 in=filt_prim_12nt_Lev_3_Tm_38_44_GC_45_55_SD_2_trim.fasta outm=cp1536.sam outu=cp1536.unaligned.sam perfectmode=t
blat seq_files/bc25mer.240k.2bit seq_files/filt_prim_12nt_Lev_3_Tm_38_44_GC_45_55_SD_2_trim.fasta BCs_psl/bc1536.psl -tileSize=12 -stepSize=1 -minMatch=1 -minScore=0 -minIdentity=100 -noHead -t=dna -q=dna

sleep 8

pslPretty BCs_psl/bc1536.psl seq_files/bc25mer.240k.2bit seq_files/filt_prim_12nt_Lev_3_Tm_38_44_GC_45_55_SD_2_trim.fasta BCs_psl/bc1536_pretty.txt

