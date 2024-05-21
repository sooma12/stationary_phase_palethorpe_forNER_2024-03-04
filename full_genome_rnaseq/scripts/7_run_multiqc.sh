#!/bin/bash

source ./config.cfg

module load anaconda3/2021.05
eval "$(conda shell.bash hook)"
conda activate /work/geisingerlab/conda_env/multiQC

multiqc $BASE_DIR --ignore data/ --ignore software/ --ignore ref/
