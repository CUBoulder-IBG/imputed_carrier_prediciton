#!/bin/bash
#SBATCH --partition=ssky-preemptable
#SBATCH --mem=2gb
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=1-13
#SBATCH --cpu-freq=2500
#SBATCH -J ROC
#SBATCH -o /work/KellerLab/meng/acc/ROC.out
#SBATCH -e /work/KellerLab/meng/acc/ROC.err

sub=$SLURM_ARRAY_TASK_ID
n1=$(expr \( $sub - 1 \) \* 10 + 1)
n2=$(expr $sub \* 10)

ml R
if [ $sub -gt 12 ]; then
Rscript /work/KellerLab/meng/acc/ROC.r 121 123
else
Rscript /work/KellerLab/meng/acc/ROC.r "$n1" "$n2"
fi





