#!/usr/bin/env bash

#! RUN : bash generatePanel_wrapper.sh <TAG> <CHR>

CHR=$1
PANELS=$2

[[ -z "$CHR" ]] && { echo "ERROR: No CHR provided for this run"; exit 1; }
[[ -z "$PANELS" ]] && { echo "ERROR: No PANELS provided for this run"; exit 1; }

JOBS=""

for TAG in `echo $PANELS | tr ':' "\n"`; do
  source ${TAG}.config

  mkdir -p "$TAG/logs"; cd $TAG
  cp ../${TAG}.config .

  VCF=`dirname $VCF_ALL`
  VCF+="broad-chr${CHR}.vcf.gz"

  # prepare files
  # Take VCF value and look to see if need to make CHR specific VCF file
  if [[ ! -e $VCF_ALL && -e $VCF ]]
  then 
    echo "ERROR - Unable to find VCF file for ${TAG}";
    exit 1;
  fi

  # The construct ${jid1##* } isolates the last word from "Submitted batch job XXXXXX" - returns XXXXXX
  # vcfChr.sh
  jid1=$(sbatch -J ${TAG}.vcfChr ${HOME}/scripts/Imputation/slurm/vcfChr.sh ${TAG} ${CHR})
  # vcfQC.sh - convert VCF to plink; perform QC; check SNPs; convert back to VCF
  jid2=$(sbatch -J ${TAG}.vcfQC --dependency=afterok:${jid1##* } ${HOME}/scripts/Imputation/slurm/vcfQC.sh ${TAG} ${CHR})

  echo $jid1
  echo $jid2

  JOBS+="${jid2##* }:"

done

echo "${JOBS::-1}"
