#!/bin/bash

# By Ryan Abramowitz

# bsub -q long -R "rusage[mem=40000] span[hosts=1] select[rh=8]" -W 8:00 -n 8 < metadata_getter.sh



module load java/1.8.0_171
java -Xmx2000m -jar /project/umw_rene_maehr/programs/Drop-seq_tools-2.5.1/3rdParty/picard/picard.jar NormalizeFasta INPUT=Homo_sapiens.GRCh38.dna.primary_assembly.fa OUTPUT=./dropseq_metadata.fasta
java -Xmx2000m -jar /project/umw_rene_maehr/programs/Drop-seq_tools-2.5.1/3rdParty/picard/picard.jar CreateSequenceDictionary REFERENCE=./dropseq_metadata.fasta OUTPUT=./dropseq_metadata.dict SPECIES=human
java -Xmx2000m -jar /project/umw_rene_maehr/programs/Drop-seq_tools-2.5.1/jar/dropseq.jar FilterGtf GTF=Homo_sapiens.GRCh38.106.filtered.gtf SEQUENCE_DICTIONARY=./dropseq_metadata.dict OUTPUT=./dropseq_metadata.gtf
java -Xmx10000m -jar /project/umw_rene_maehr/programs/Drop-seq_tools-2.5.1/jar/dropseq.jar ConvertToRefFlat ANNOTATIONS_FILE=./dropseq_metadata.gtf SEQUENCE_DICTIONARY=./dropseq_metadata.dict OUTPUT=./dropseq_metadata.refFlat
java -Xmx10000m -jar /project/umw_rene_maehr/programs/Drop-seq_tools-2.5.1/jar/dropseq.jar ReduceGtf GTF=./dropseq_metadata.gtf SEQUENCE_DICTIONARY=./dropseq_metadata.dict OUTPUT=./dropseq_metadata.reduced.gtf
java -Xmx10000m -jar /project/umw_rene_maehr/programs/Drop-seq_tools-2.5.1/jar/dropseq.jar CreateIntervalsFiles SEQUENCE_DICTIONARY=./dropseq_metadata.dict REDUCED_GTF=./dropseq_metadata.reduced.gtf PREFIX=dropseq_metadata OUTPUT=.
#mkdir -p ./STAR
module load star/2.7.10a
/share/pkg/star/2.7.10a/bin/STAR --runMode genomeGenerate --genomeDir ./STAR --genomeFastaFiles ./dropseq_metadata.fasta --sjdbGTFfile ./dropseq_metadata.gtf --sjdbOverhang 100 --runThreadN 8
module load samtools/1.15.1
bgzip ./dropseq_metadata.fasta
samtools faidx ./dropseq_metadata.fasta.gz

mv ./STAR ./STAR_index
gunzip dropseq_metadata.fasta.gz

