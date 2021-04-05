#!/usr/bin/env bash

#! RUN : sbatch vcfQC.sh <TAG> <CHR>

#! sbatch directives begin here ###############################
#! Which project should be charged:
#SBATCH -A GODOGS-SL2-CPU
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 1:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o logs/job-%j.out

module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

module load bcftools-1.9-gcc-5.4.0-b2hdt5n        # bcftools
module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7    # bgzip/tabix

PROJECT=$1
source setup.config

REF=""
for p in `echo $REF_PANELS | tr ':' "\n"`; do
  REF+="$p/$p.vcf.gz "
done
echo $REF

bcftools merge -m id --force-samples --missing-to-ref ${REF} > merged.vcf
bcftools query -l merged.vcf | grep ':' > samples.dups
bcftools view -Ov -S ^samples.dups merged.vcf | bcftools view -m2 -M2 -v snps | bgzip -c > ${PROJECT}.ref-panel.vcf.gz
tabix -p vcf ${PROJECT}.ref-panel.vcf.gz

rm -rf merged.vcf samples.dups