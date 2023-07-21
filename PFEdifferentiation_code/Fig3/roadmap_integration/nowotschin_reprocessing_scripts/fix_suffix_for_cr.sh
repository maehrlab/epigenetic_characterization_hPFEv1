# symlink fq's into desired naming format: "MySample_S1_L001_I1_001.fastq.gz"
#
for folder in $(ls . | grep -e 'E.\.[75|5]' ); do
  for read_num in 1 2; do
    ii=0
    for fq in ${folder}/*_${read_num}.fastq.gz; do
      fq_basename=$(basename $fq)
      fq_dirname=$(dirname $fq)
      #echo $fq_basename
      #echo $fq_dirname
      #echo $SRR_num
      lane=$(($(($ii / 4)) + 1))
      index=$(($(($ii % 4)) + 1))
      ln -bs ${fq_basename} ${fq_dirname}/SRR_S1_L00${lane}_R${read_num}_00${index}.fastq.gz
      ii=$((ii + 1)) 
    done 
  done
done
