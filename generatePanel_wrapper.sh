#!/usr/bin/env bash

#! RUN : bash generatePanel_wrapper.sh <CHR> <PROJECT>
#! Eg. : bash generatePanel_wrapper.sh 13 go_lab

CHR=$1
PROJECT=$2
SCRIPTS=`dirname $0`

# Check that variables are provided and if not exit the script
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
  jid1=$(sbatch -J ${TAG}.vcfChr ${SCRIPTS}/slurm/vcfChr.sh ${TAG} ${CHR})
  # vcfQC.sh - convert VCF to plink; perform QC; check SNPs; convert back to VCF
  jid2=$(sbatch -J ${TAG}.vcfQC --dependency=afterok:${jid1##* } ${SCRIPTS}/slurm/vcfQC.sh ${TAG} ${CHR})

  #echo $jid1
  #echo $jid2

  JOBS+="${jid2##* }:"

  cd ../

done

# prepare BAM files...
if [ ! -e $BAM_FILES/chr${CHR} ]
then
  echo "Need to create chromosome speciic BAM files - $BAM_FILES";
  mkdir "${BAM_FILES}/chr${CHR}"
  # SUBMIT JOB ARRAY!
  count=`ls ${BAM_FILES}/*.bam | wc -l`
  jid3=$(sbatch -J ${PROJECT}.bams --array=1-${count} ${SCRIPTS}/slurm/prepareBamFiles.sh ${PROJECT} ${CHR})
  JOBS+="${jid3##* }:"
fi


# once all 2nd jobs AND job#3 have finished
# merge reference panles together - details given by REF_PANELS varaible in setup.config

# mergeRefPanels.sh
jid4=$(sbatch -J ${PROJECT}.mergeRef --dependency=afterok:${JOBS::-1} ${SCRIPTS}/slurm/mergeRefPanels.sh ${PROJECT} ${CHR})

# fillProjectGaps.sh
jid5=$(sbatch -J ${PROJECT}.fillGaps --array=1-16 --dependency=afterok:${jid4##* } ${SCRIPTS}/slurm/fillProjectGaps.sh ${PROJECT} ${CHR})

# addRefSNPs.sh
jid6=$(sbatch -J ${PROJECT}.addRefSNPs --dependency=afterok:${jid5##* } ${SCRIPTS}/slurm/addRefSNPs.sh ${PROJECT} ${CHR})

# preparePanel.sh
jid7=$(sbatch -J ${PROJECT}.preparePanel --dependency=afterok:${jid6##* } ${SCRIPTS}/slurm/preparePanel.sh ${PROJECT} ${CHR})

# shapeitPanel.sh
sbatch -J ${PROJECT}.shapeit --dependency=afterok:${jid7##* } ${SCRIPTS}/slurm/shapeitPanel.sh ${PROJECT} ${CHR}
