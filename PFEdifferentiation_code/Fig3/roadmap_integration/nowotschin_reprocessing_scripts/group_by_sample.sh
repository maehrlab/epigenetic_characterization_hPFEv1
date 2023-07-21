dos2unix ../metadata/*
for sample in $(cut -f1 -d',' ../metadata/peer_fq_by_sample.csv | tail -n +2)
do
  mkdir -p $sample
done

for line in $(cat ../metadata/peer_sample_by_fq.csv | tail -n +2 )
do 
  sample=$(echo $line | cut -f2 -d',' )
   fastq=$(echo $line | cut -f1 -d',' )
  mv ${fastq}_1.fastq.gz  $sample
  mv ${fastq}_2.fastq.gz  $sample
done