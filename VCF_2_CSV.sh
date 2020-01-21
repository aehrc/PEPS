#!/bin/bash

VcfFile=$1          # Pass to the VCF file

echo "=============================================================="

echo " Argument 1: path to vcf or vcf.bgz file"

echo "=============================================================="

echo " Argument 1: $VcfFile"

echo "=============================================================="

echo " ***** This scripts uses bcftools."
echo " ***** This scripts uses vcf file CHROM:POS:REF:ALT as SNP id."
echo " ***** The output is vcffile.csv"
echo " ***** This scripts does not check input arguments."


echo " >>>>> Converting...."

bcftools query -f '%CHROM:%POS:%REF:%ALT\n' $VcfFile > $VcfFile.id.tsv
bcftools query -l  $VcfFile | awk 'BEGIN{print("SNP")}{print}' | tr \\n , | sed 's/.$//g' | awk '{print}' > $VcfFile.head
bcftools query -f '[%GT\t]\n' $VcfFile | tr -d '|' | tr -d '/' | tr '\.' '0' | sed 's/00/0/g' | sed 's/01/1/g' | sed 's/10/1/g' | sed 's/11/2/g' | paste $VcfFile.id.tsv - | tr \\t ',' | sed 's/.$//g' | cat $VcfFile.head - > $VcfFile.csv

echo " >>>>> Done"

echo " >>>>> Remove temporary files"
rm $VcfFile.id.tsv $VcfFile.head

echo " >>>>> Output: $VcfFile.csv"
