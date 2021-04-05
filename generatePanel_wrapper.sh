#!/usr/bin/env bash

#! RUN : bash generatePanel_wrapper.sh <CHR> <PROJECT>
#! Eg. : bash generatePanel_wrapper.sh 13 go_lab

CHR=$1
PROJECT=$2

[[ -z "$CHR" ]] && { echo "ERROR: No CHR provided for this run"; exit 1; }
[[ -z "$PROJECT" ]] && { echo "ERROR: No PROJECT provided for this run"; exit 1; }

JOBS=""
mkdir -p chr${CHR}/logs; cd chr${CHR};
cp ../setup.config .
source setup.config

for TAG in `echo $REF_PANELS:$PROJECT | tr ':' "\n"`; do

  echo $TAG

  mkdir "$TAG"; cd $TAG
  cp ../setup.config .

  VCF_ALL=${!TAG}
  VCF=`dirname $VCF_ALL`
  VCF+="${TAG}-chr${CHR}.vcf.gz"

  # prepare files
  # Take VCF value and look to see if need to make CHR specific VCF file
  if [[ ! -e $VCF_ALL && -e $VCF ]]
  then 
    echo "ERROR - Unable to find VCF file for ${TAG}";
    exit 1;
  fi

  # The construct ${jid1##* } isolates the last word from "Submitted batch job XXXXXX" - returns XXXXXX
  # vcfChr.sh
  #jid1=$(sbatch -J ${TAG}.vcfChr ${HOME}/scripts/Imputation/slurm/vcfChr.sh ${TAG} ${CHR})
  # vcfQC.sh - convert VCF to plink; perform QC; check SNPs; convert back to VCF
  #jid2=$(sbatch -J ${TAG}.vcfQC --dependency=afterok:${jid1##* } ${HOME}/scripts/Imputation/slurm/vcfQC.sh ${TAG} ${CHR})

  #echo $jid1
  #echo $jid2

  #JOBS+="${jid2##* }:"

  cd ../

done

#echo "${JOBS::-1}"

pwd

# once all 2nd jobs have finished
# merge reference panles together - details by REF_PANELS varaible in setup.config
# mergeRefPanels.sh
# jid3=$(sbatch -J ${PROJECT}.mergeRed --dependency=afterok:${JOBS::-1} ${HOME}/scripts/Imputation/slurm/mergeRefPanels.sh ${PROJECT})
jid3=$(sbatch -J ${PROJECT}.mergeRef ${HOME}/scripts/Imputation/slurm/mergeRefPanels.sh ${PROJECT})
echo $jid3
