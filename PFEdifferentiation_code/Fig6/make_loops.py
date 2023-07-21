#!/usr/bin/env python

# Code by: Ryan Abramowitz
from celloracle import motif_analysis as ma
import os

# conda activate celloracle

basedir = './'
for sample in ['hESC','hDE','hAFE','hPFE']:
    bed = ma.read_bed(basedir + sample + f'/{sample}_activator_regions_with_footprints_tmp.bed')
    peaks = ma.process_bed_file.df_to_list_peakstr(bed)
    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38") # comes with cell_oricle; if peak contains tss, retain peek and add gene. Else, delete peak
    tss_annotated.to_csv(basedir + sample +"/tss_peaks_for_looping.bed",sep='\t',index=None,header=None)
    os.system(f"""bedtools window -a {basedir+sample}/{sample}_activator_regions_with_footprints_tmp.bed -b {basedir+sample}/tss_peaks_for_looping.bed -w 100000 > {basedir+sample}/{sample}_loops.bed""")
