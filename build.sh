#!/bin/bash

# Run the pipeline on the cluster
#srun -c2 --mem 80G -J build_haemdata-Rpkg 
/opt/singularity/3.7.0/bin/singularity exec \
    --cleanenv -H $HOME/.singularity/home:/home/$USER \
    /opt/singularity-images/kmorris/rbioc_3.14.006.sif \
    R CMD BATCH build.R  > /dev/null 2>&1 &

# Alternately, run the pipeline in a persistent background process:
# nohup nice -4 R CMD BATCH run.R &

# Removing .RData is recommended.
rm -f .RData
