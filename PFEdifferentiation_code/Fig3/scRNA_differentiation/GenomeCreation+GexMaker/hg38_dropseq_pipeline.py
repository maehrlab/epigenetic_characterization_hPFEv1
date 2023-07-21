#!/bin/env python

# By Ryan Abramowitz

from __future__ import print_function
import argparse
import os
import subprocess
import datetime
# scp -i ~/sshkey-rsa4096 dropseq_generate_dge.py ra73w@ghpcc06.umassrc.org:/project/umw_rene_maehr/programs/hg38_dropseq_pipeline.py
# ssh ra73w@ghpcc06.umassrc.org -i ~/sshkey-rsa4096
# bsub -q interactive -R "rusage[mem=40000] span[hosts=1] select[rh=8]" -Is -W 8:00 -n 8 bash
# cd /nl/umw_rene_maehr/RA_projects/redo
# /project/umw_rene_maehr/programs/hg38_dropseq_pipeline.py --fastq1 /nl/umw_rene_maehr/backup/process/2015DEC10_PE_20_50_Dropseq/hAFE_H1_rep1_2015DEC10Dropseq.1.fastq.gz --fastq2 /nl/umw_rene_maehr/backup/process/2015DEC10_PE_20_50_Dropseq/hAFE_H1_rep1_2015DEC10Dropseq.2.fastq.gz --sample_id sample0 --output .
# find . -name "*.bam" -type f -delete
# find . -name "*" -type f | grep 'dge\|readcount'
# find . -name "*" -type f | grep 'dge\|readcount' | xargs -I file echo $PWD'/file'
# find . -name "*" -type f | grep -v 'dge\|readcount' | xargs -I file rm file

parser = argparse.ArgumentParser()
parser.add_argument('--fastq1', help = 'FASTQ file with barcode reads. Specify this or --tagged_bam.')
parser.add_argument('--fastq2', help = 'FASTQ file with exon reads.    Specify this or --tagged_bam.')
parser.add_argument('--tagged_bam', help = 'Barcode-tagged, unaligned BAM file to begin with. Specify this or the FASTQ inputs.')
parser.add_argument('--sample_id', help = 'Used as a results folder name and a file prefix.')
parser.add_argument('--output', help = 'Results go hear with additional layers for date and sample id.')
parser.add_argument('--test_mode', help = 'y/n. Subset a thousand reads for a quick test run. Appends _TEST to sample id.',
                    default = 'n' )
parser.add_argument('--min_genes_per_cell', help = 'Passed to DigitalGeneExpression.',
                    default = 1000, type = int)
parser.add_argument("--n_thread", help="Goes into bsub's -n. Usually we use 10.",
                    default=10, type=int)
parser.add_argument("--runtime_cap", help="Goes into bsub's -W. Usually we use 10:00",
                    default="10:00")
parser.add_argument("--mem", help="Goes into bsub's -R rusage[mem=<mem>]. Usually we use 32000 because STAR is a glutton for RAM.",
                    default=32000, type=int)
parser.add_argument("--space_saver", help="y/n. If yes, remove the largest intermediate files.",
                    default="y")
parser.add_argument("--queue", help="LSF queue to use.",
                    default="long")

args = parser.parse_args()

############ Parse inputs

# Test mode
if args.test_mode == 'y':
    print ('Running in test mode.')
    args.sample_id += '_TEST'
    # if args.tagged_bam is None:
    #     raise Exception('Sorry, test mode is currently available only for tagged BAM input.')

# Reference genome
GENOME_DIR = '/project/umw_rene_maehr/genomes/Homo_sapiens/hg38/dropseq'
REFFLAT    = '/project/umw_rene_maehr/genomes/Homo_sapiens/hg38/dropseq/dropseq_metadata.refFlat'

############ Set up output folder and filenames
TODAY = datetime.datetime.now().strftime('%Y_%m_%d__%H_%M_%S')
OUTPUT = os.path.join( args.output, args.sample_id, TODAY )
os.makedirs( OUTPUT, mode = 0o744 )

non_star_prefix = os.path.join( OUTPUT, args.sample_id )
SUMMARY_TAG_AAAA = non_star_prefix + '_polyA_summary.txt'
SUMMARY_ADAPTER  = non_star_prefix + '_adapter_trimming_report.txt'
SUMMARY_TAG_CELL = non_star_prefix + '_tagged_cell_summary.txt'
SUMMARY_TAG_UMIS = non_star_prefix + '_tagged_UMIs_summary.txt'
FASTQ_CLEAN      = non_star_prefix + '_tagged_filtered_polyA.fastq'
TEST_BAM         = non_star_prefix + '_shortened_test_input.bam'

STAR_PREFIX = os.path.join( OUTPUT, 'star_' + args.sample_id )
STAR_SAM_ORIG = STAR_PREFIX + 'Aligned.out.sam'
STAR_SAM_SORT = STAR_PREFIX + 'Aligned_sorted.bam'
STAR_EXON_TAG = STAR_PREFIX + '_merged_exon_tagged.bam'
STAR_EXON_CLE = STAR_PREFIX + '_merged_exon_tagged_clean.bam'
STAR_EXON_CLR = STAR_PREFIX + '_merged_exon_tagged_cleaner.bam'
STAR_BEAD_SYN = STAR_PREFIX + '_bead_synthesis_stats.txt'
STAR_BEAD_SUM = STAR_PREFIX + '_bead_synthesis_stats_summary.txt'
STAR_READ_CNT = STAR_PREFIX + '_merged_gene_exon_tagged_clean_readcount.txt.gz'
STAR_EXON_DGE = STAR_PREFIX + '_merged_gene_exon_tagged_clean.dge.txt.gz'
DGE_SUMMARY   = STAR_PREFIX + '_merged_gene_exon_tagged_clean.dge.summary.txt'

TEMPDIR = os.path.join( OUTPUT, args.sample_id, 'temp' )
os.makedirs( TEMPDIR, mode = 0o744 )
tempbam = lambda x: os.path.join( TEMPDIR,'temp' + str(x) + '.bam')

# Option to skip cell/UMI tagging by starting with a tagged BAM
if args.tagged_bam is not None:
    BAM_TAG_AAAA = args.tagged_bam
    print('Beginning with tagged BAM rather than paired FASTQ files.')
    if ( args.fastq1 is not None ) or ( args.fastq2 is not None ):
        raise Exception('You cannot give FASTQ input and also a tagged bam.' + \
                        'Solution: Specify --tagged_bam or --fastq1 and --fastq2 but not both.')
else:
    BAM_TAG_AAAA = os.path.join( OUTPUT, args.sample_id + '_tagged_filtered_polyA.bam' )

print('Results will be sent to ' + OUTPUT)
print('Genome used: ' + GENOME_DIR)
print('Exons used:  ' + REFFLAT)
print('Min genes per cell: %d%', args.min_genes_per_cell)
print ('\n'.join(['================ Files to expect: ==================\n',
                 SUMMARY_TAG_AAAA,
                 SUMMARY_TAG_CELL,
                 SUMMARY_ADAPTER,
                 SUMMARY_TAG_UMIS,
                 BAM_TAG_AAAA,
                 FASTQ_CLEAN,
                 STAR_SAM_ORIG,
                 STAR_SAM_SORT,
                 STAR_EXON_TAG,
                 STAR_EXON_CLE,
                 STAR_EXON_CLR,
                 STAR_BEAD_SYN,
                 STAR_BEAD_SUM,
                 STAR_READ_CNT,
                 STAR_EXON_DGE,
                 DGE_SUMMARY,
                 '\n================                  ==================']) )

############ Set up shorthand for commands
# There's no input/output in this section.
# To keep the key files separate from the clutter of options shown here, the input/output is shown below.

# PICARD = 'java -Xmx6g -jar /share/pkg/picard/1.96/'
# DROPSEQ_TOOLS = '/project/umw_rene_maehr/programs/Drop-seq_tools-1.0/'
# STAR = '/share/pkg/star/2.4.2a/STAR'
PICARD = 'java -Xmx6g -jar /project/umw_rene_maehr/programs/Drop-seq_tools-2.5.1/3rdParty/picard/picard.jar '
DROPSEQ_TOOLS = '/project/umw_rene_maehr/programs/Drop-seq_tools-2.5.1/'
STAR = '/share/pkg/star/2.7.10a/bin/STAR'

# Convert fastq to sam
fastq_to_bam = PICARD + 'FastqToSam SM=' + args.sample_id
# tag unaligned BAM file with Cell Barcode
barcode_cell = DROPSEQ_TOOLS + 'TagBamWithReadSequenceExtended SUMMARY=' + SUMMARY_TAG_CELL + \
               ' BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1'
# tag unaligned BAM file with Unique Molecular Identifier (UMI)
barcode_mol = DROPSEQ_TOOLS + 'TagBamWithReadSequenceExtended SUMMARY=' + SUMMARY_TAG_CELL + \
              ' BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1'
# remove reads with low quality cellular barcode
filter_cells = DROPSEQ_TOOLS + 'FilterBam TAG_REJECT=XQ'
# trim the SMART Adapter that can occur 5' of the read
trim_start = DROPSEQ_TOOLS + 'TrimStartingSequence OUTPUT_SUMMARY=' + SUMMARY_ADAPTER + \
             ' SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5'
# polyA trimmer
poly_a_trim = DROPSEQ_TOOLS + 'PolyATrimmer OUTPUT_SUMMARY=' + SUMMARY_TAG_AAAA + ' MISMATCHES=0 NUM_BASES=6'
# convert cleaned, tagged and unaligned BAM to FASTQ file with Picard tools
bam_to_fastq = PICARD + 'SamToFastq'
# STAR alignment (eventually with random choice of primary read from among multiple equally good alignments, but not working yet)
run_star = STAR + ' --genomeDir ' + os.path.join( GENOME_DIR, 'STAR_index') + ' --runThreadN 10 '  # --outMultimapperOrder Random --outSAMmultNmax 1'
# bam sorting for Star alignment results
sort_star_sam = PICARD + 'SortSam SO=queryname'
# merge bam alignment from star
merge_bams =    PICARD + 'MergeBamAlignment REFERENCE_SEQUENCE=' + os.path.join( GENOME_DIR, 'dropseq_metadata.fasta' ) + \
                '  INCLUDE_SECONDARY_ALIGNMENTS=False PAIRED_RUN=False'
# tag reads with exons
tag_exons = DROPSEQ_TOOLS + 'TagReadWithGeneFunction ANNOTATIONS_FILE=' + REFFLAT
# DetectBeadSynthesisError
bead_check = DROPSEQ_TOOLS + \
             'DetectBeadSynthesisErrors OUTPUT_STATS=' + STAR_BEAD_SYN + \
             ' SUMMARY=' + STAR_BEAD_SUM + \
             ' NUM_BARCODES=2000 PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC TMP_DIR=' + TEMPDIR
# DetectBeadSubstitutionErrors
bead_check2 = DROPSEQ_TOOLS + 'DetectBeadSubstitutionErrors -TMP_DIR' + TEMPDIR


# Get cellwise readcounts
read_counts = DROPSEQ_TOOLS + 'BamTagHistogram TAG=XC'
# Digital gene expression for cells with >1000 genes
make_dge = DROPSEQ_TOOLS + 'DigitalExpression SUMMARY=' + DGE_SUMMARY + ' TMP_DIR=' + TEMPDIR + \
    ' MIN_NUM_GENES_PER_CELL=' + str(args.min_genes_per_cell)

############ Form pipeline, showing inputs and outputs.

# Do cell/molecular barcoding unless starting with a tagged BAM
if args.tagged_bam is None:
    barcode = fastq_to_bam + \
              ' F1=' + str(args.fastq1) + \
              ' F2=' + str(args.fastq2) + \
              ' OUTPUT='+tempbam(0)+'\n\n' + \
              barcode_cell + ' INPUT='+tempbam(0)+' OUTPUT='+tempbam(1)+'\n\n' + \
              barcode_mol + ' INPUT='+tempbam(1)+' OUTPUT='+tempbam(2)+'\n\n' + \
              filter_cells + ' INPUT='+tempbam(2)+' OUTPUT='+tempbam(3)+'\n\n' + \
              trim_start + ' INPUT='+tempbam(3)+' OUTPUT='+tempbam(4)+'\n\n' + \
              poly_a_trim + ' INPUT='+tempbam(4)+' OUTPUT=' + BAM_TAG_AAAA
else:
    barcode = ""

# bam to fastq with option for small test run
if args.test_mode == 'y':
    truncate_bam = 'samtools view -h ' + BAM_TAG_AAAA + ' | head -n 1000 | samtools view -bS > ' + TEST_BAM
    bam2fastq = truncate_bam + '\n\n' + bam_to_fastq + ' INPUT=' + TEST_BAM + ' FASTQ=' + FASTQ_CLEAN
else:
    bam2fastq = bam_to_fastq + ' INPUT=' + BAM_TAG_AAAA + ' FASTQ=' + FASTQ_CLEAN

# Alignment
run_star += ' --readFilesIn ' + FASTQ_CLEAN + ' --outFileNamePrefix ' + STAR_PREFIX

# sorting
sort_star_sam+=' I=' + STAR_SAM_ORIG + ' O=' + STAR_SAM_SORT

# exon tagging
merge_and_tag = merge_bams + \
                ' UNMAPPED_BAM=' + BAM_TAG_AAAA + \
                ' ALIGNED_BAM=' + STAR_SAM_SORT + \
                ' OUTPUT='+tempbam(5)+'\n\n' + \
                tag_exons + \
                ' I='+tempbam(5)+' O=' + STAR_EXON_TAG

# Bead synthesis error checking
check_beads = bead_check + ' I=' + STAR_EXON_TAG + ' O=' + STAR_EXON_CLE
check_beads2 = bead_check + ' I=' + STAR_EXON_CLE + ' O=' + STAR_EXON_CLR
# final step
read_counts +=' I=' + STAR_EXON_CLR + ' O=' + STAR_READ_CNT
make_dge    +=' I=' + STAR_EXON_CLR + ' O=' + STAR_EXON_DGE

############ List and run commands

commands = ["#!/bin/sh",
            "module load java/1.8.0_171",
            "module load star/2.7.10a",
            "module load samtools/1.15.1",
            barcode,
            bam2fastq,
            run_star,
            sort_star_sam,
            merge_and_tag,
            check_beads,
            check_beads2,
            read_counts,
            make_dge]

# Delete largest intermediates: the uncompressed SAM from STAR and the FASTQ fed into STAR.
if args.space_saver == 'y':
    commands = commands +  ["rm " + STAR_SAM_ORIG]
    commands = commands +  ["rm " + FASTQ_CLEAN]
    commands = commands +  ["rm -r " + TEMPDIR]
    commands = commands +  ["ls "+os.path.join( OUTPUT, '*.bam' )+" | xargs rm"]

# Save and execute commands
with open(os.path.join(OUTPUT, 'commands.sh'), 'a+') as command_file:
    command_file.writelines( "\n\n".join( commands ) )
os.chmod(os.path.join(OUTPUT, 'commands.sh'), 0o755)
bsub_wrapper = ' bsub -W ' +        str(args.runtime_cap) + \
                ' -R rusage[mem=' + str(args.mem / args.n_thread) + '] '  + \
                ' -q ' +            str(args.queue) + ' '  + \
                ' -n ' +            str(args.n_thread)    + \
                ' -R span[hosts=1] '                      + \
                ' " ' + os.path.join(OUTPUT, 'commands.sh') + ' " '
with open(os.path.join(OUTPUT, 'bsub_commands.sh'), 'a+') as bsub_file:
    bsub_file.writelines( "\n".join(["#!/bin/sh", bsub_wrapper]) )
os.chmod(os.path.join(OUTPUT, 'bsub_commands.sh'), 0o755)
subprocess.call( os.path.join(OUTPUT, 'bsub_commands.sh') )


