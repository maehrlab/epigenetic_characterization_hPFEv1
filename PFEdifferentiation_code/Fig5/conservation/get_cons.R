# Code By: Jack Huey
library(rGREAT)
library(rtracklayer)
library(parallel)

ont=c("GO Molecular Function","GO Biological Process","GO Cellular Component","Human Phenotype","Ensembl Genes")

do.ontology = function(filepath) {
  tmp = import.bed(filepath)
  if (length(tmp) == 0) {
    return()
  }
  if (file.exists(paste0(gsub(".bed", "", filepath), ".", gsub(" ", "_", ont[1]), ".csv"))) {
    return()
  }
  job = submitGreatJob(tmp,species='hg38')
  tb1 = getEnrichmentTables(job,ontology = ont)
  for (info in ont){
    write.csv(data.frame(tb1[info]), paste0(gsub(".bed", "", filepath), ".", gsub(" ", "_", info), ".csv"))
  }
}

mclapply(list.files("output/tmp", pattern = ".conserved_peaks.bed", full.names = T), do.ontology, mc.cores = 2)
