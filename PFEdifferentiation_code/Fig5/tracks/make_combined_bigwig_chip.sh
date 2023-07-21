#!/bin/bash
. "/home/ml33w/miniconda3/etc/profile.d/conda.sh"

#module load samtools/1.9
#conda init bash
conda activate /home/ml33w/miniconda3/envs/deeptools


basedir="/nl/umw_rene_maehr/JH_projects/chip/grn-1/"
samp="${1}"
outdir="/nl/umw_rene_maehr/ML_projects/bulkchip_processing/"


bw1=${basedir}/${samp}/out/signal/${2}/$(ls ""${basedir}/${samp}/out/signal/${2}/"" | grep "fc")
bw2=${basedir}/${samp}/out/signal/${3}/$(ls ""${basedir}/${samp}/out/signal/${3}/"" | grep "fc")

op="${outdir}/${samp}_merged_fc"



bigwigCompare -b1 ${bw1} -b2 ${bw2} --skipZeroOverZero --operation mean  --outFileFormat bigwig --numberOfProcessors max -o ${op}.bw