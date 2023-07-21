

# Code by: Ryan Abramowitz
cd /home/ryan/Desktop/GRN/GRN_Ryan_ATAC_only/

for VARIABLE in hPFE hAFE hDE hESC
do

mkdir ${VARIABLE}

bedtools intersect -wa -wb -a ../ATAC/${VARIABLE}.bed -b ../ATAC/final_peakset_fp.bed > ${VARIABLE}/${VARIABLE}_activator_regions_with_footprints.bed

cat ${VARIABLE}/${VARIABLE}_activator_regions_with_footprints.bed | awk '{print $1"\t"$2"\t"$3}' | bedtools sort | uniq > ${VARIABLE}/${VARIABLE}_activator_regions_with_footprints_tmp.bed

done





