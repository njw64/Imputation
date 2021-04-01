#!/usr/bin/env bash

#! RUN : bash generatePanel_wrapper.sh <TAG> <CHR>

TAG=$1
CHR=$2

[[ -z "$TAG" ]] && { echo "ERROR: No TAG given for this run"; exit 1; }
[[ -z "$CHR" ]] && { echo "ERROR: No CHR provided for this run"; exit 1; }

source ${TAG}.config
module load bcftools-1.9-gcc-5.4.0-b2hdt5n
module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7

# prepare files
# Take VCF value and look to see if need to make CHR specific VCF file

if [[ ! -e $VCF_ALL && -e $VCF ]]
then 
  echo "ERROR - Unable to find VCF file for ${TAG}";
  exit 1;
fi


# slurm_vcfChr.sh
tabix -p vcf $VCF_ALL
if [ ! -e $VCF ]
then
  echo "Need to create $VCF";
  echo "chr$CHR $CHR" > $TAG.chrs
  bcftools view $VCF_ALL --regions $CHR,chr$CHR | bcftools view -m2 -M2 -v snps | bcftools filter -e "QUAL < 20" | bcftools annotate --rename-chrs $TAG.chrs | bcftools annotate --set-id +'%CHROM:%POS' | bgzip -c > $VCF
  tabix -p vcf $VCF
fi

echo $VCF
echo $VCF_ALL