#!/usr/bin/env bash

#! RUN : sbatch prepareGWAS.sh <TAG> <CHR>

#! sbatch directives begin here ###############################
#! Which project should be charged:
#SBATCH -A GODOGS-SL2-CPU
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 0:30:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o ../logs/job-%j.out

module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

module load bcftools-1.9-gcc-5.4.0-b2hdt5n        # bcftools
module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7    # bgzip/tabix
module load plink-1.9-gcc-5.4.0-sm3ojoi           # plink/plinkdog

CHR=$1
source ../setup.config

# QC of gwas data - mid/geno/maf
# Pull out the correct chromosome only
CMD="${PLINK} --chr ${CHR} --bfile ${GWAS_DIR}/${GWAS_PLINK} --snps-only --make-bed"

if [ ! -z "$SAMPLES" ]
then
  if [ ! -e $SAMPLES ]
  then
    echo "ERROR - Samples file provided but not found [${SAMPLES}]";
    #exit 1;
  fi
  
  # Filter GWAS data by samples as well as CHR
  CMD+=" --keep ${SAMPLES}"
fi

eval $CMD

perl -lane 'if(($F[4] eq "A" && $F[5] eq "T") || ($F[4] eq "T" && $F[5] eq "A") || ($F[4] eq "C" && $F[5] eq "G") || ($F[4] eq "G" && $F[5] eq "C")){ print $F[1] }' plink.bim > ambiguous.snps
eval "$PLINK --bfile plink --geno 0.03 --mind 0.05 --maf 0.01 --hwe 0.00005 --freq --missing --make-bed --out ${GWAS_PLINK}"
