    mm10_cr=/project/umw_rene_maehr/genomes/refdata-cellranger-mm10-3.0.0
cellranger3=/project/umw_rene_maehr/programs/cellranger-3.0.1/cellranger
for sample in $(cut -f1 -d, ../metadata/peer_fq_by_sample.csv) 
do  
   echo bsub -W 20:00 -q long "${cellranger3} count \
                   --id=$(echo $sample | sed s/\\./_/ ) \
                   --sample=SRR \
                   --transcriptome=${mm10_cr} \
                   --fastqs=../fastq/${sample} \
                   --chemistry=SC3Pv2 \
                   --jobmode=lsf"
done
