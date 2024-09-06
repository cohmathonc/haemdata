#!/bin/bash
#SBATCH --job-name=targets
#SBATCH -c 2
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

/opt/singularity/3.7.0/bin/singularity exec \
    --env OMP_NUM_THREADS=1 \
    --env R_LIBS_SITE=/opt/singularity-images/rbioc/rlibs/bioc-3.18 \
    -B /opt,/labs,/run \
    /opt/singularity-images/rbioc/vscode-rbioc_3.18.sif \
    Rscript -e "targets::tar_make()"

