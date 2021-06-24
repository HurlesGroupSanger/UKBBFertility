#!/software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript --vanilla

library(data.table)
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly = T)

variant.maf <- as.double(args[1])
variant.type <- args[2]
sex <- as.integer(args[3])
outfile <- args[4]

# Read all possible genes
genes.hg38 <- fread("rawdata/genelists/hg38.coordinates.txt",fill = T)
setnames(genes.hg38,names(genes.hg38),c("GENE","CHR","START","STOP","HGNC","NAME"))
genes.hg38 <- genes.hg38[CHR != "X" & CHR != "Y"]

## Load actual SNV data:
load.data <- function(counts.file,
                      samples.file,
                      indv.to.keep) {
  
  ## Read the counts file:
  variant.counts <- fread(counts.file)
  setnames(variant.counts,names(variant.counts),c("eid","gt","chr","pos","ref","alt","CADD","MPC","gnomadAF","VQSR","PEXT","CSQ","UKBB.AC","AN","UKBB.AF","AP","GENE","last_exon","last_intron"))
  variant.counts[,eid:=as.character(eid)]
  
  ## Set a variable in the table for SNV or INDEL
  variant.counts[,allele.type:=ifelse(nchar(ref)==nchar(alt),"SNV","INDEL")]
  
  ## Set gnomadAF to 0 if it was not found for a particular variant
  variant.counts[,gnomadAF:=ifelse(is.nan(gnomadAF),0,gnomadAF)]
  
  ## Need this because some samples might have 0 variants. (SYN should always have some...?)
  ## Will load a list of samples that we should expect to have and attaches the phenotype ID to them
  final.samples <- fread(samples.file, header = F)
  setnames(final.samples,"V1","eid")
  final.samples[,eid:=as.character(eid)]
  
  # Exclude/Include individuals/variants based on ancestry and relatedness
  variant.counts <- variant.counts[eid %in% indv.to.keep]
  final.samples <- final.samples[eid %in% indv.to.keep,c("eid")]
  
  ###!! FILTERING !!###
  ## Filter by CADD & MPC:
  variant.counts <- variant.counts[!is.na(CADD)]
  variant.counts <- variant.counts[((CSQ == "LOF_HC" | (CSQ == "MIS" & MPC > 2)) & CADD > 25) | (CSQ != "LOF_HC" & CSQ != "MIS")]
  
  ## Get rid of last exon LoFs:
  variant.counts <- variant.counts[(CSQ == "LOF_HC" & last_exon == F & last_intron == F) | (CSQ != "LOF_HC")]
  
  ## Filter by gnomAD/UKB MAF:
  variant.counts <- variant.counts[gnomadAF < 0.001 & UKBB.AF < 0.001]

  ## Filter by allele missingness:
  variant.counts <- variant.counts[AN >= AP/2]
  ###!! FILTERING !!###
  
  return(list(variant.counts,final.samples))
  
}

# UKBB/Hg38
counts.file <- "rawdata/snvresources/counts.ukbb_wes.200k.txt"
samples.file <- "rawdata/snvresources/ukbb_ids.200k.txt"
indv.to.keep <- fread("rawdata/snvresources/indv_to_keep.txt",header = F)
indv.to.keep <- indv.to.keep[,V1]

# The weird paste is because of how sampleIDs are represented in the variants table
data.UKBB.200k <- load.data(counts.file,samples.file,indv.to.keep)
variants.UKBB.200k <- data.UKBB.200k[[1]]
samples.UKBB.200k <- data.UKBB.200k[[2]]
rm(data.UKBB.200k)

# Load the phenotype matrix:
UKBB.phenotype.data <- readRDS("rawdata/phenofiles/UKBB.phenotype.rdat")

# Building genotype matrix:
if (variant.maf == 0) {
  gt.dt <- variants.UKBB.200k[UKBB.AC == 1 & CSQ == variant.type & !is.na(GENE) & chr != "X" & chr != "Y",sum(gt),by=c("GENE","eid")]  
} else {
  gt.dt <- variants.UKBB.200k[UKBB.AF < variant.maf & CSQ == variant.type & !is.na(GENE) & chr != "X" & chr != "Y",sum(gt),by=c("GENE","eid")]  
}
gt.dt[,V1:=if_else(V1>1,1L,V1)]

# Building covariate matrix:
if (sex == 1) {
  y.var <- "children.fathered"
} else {
  y.var <- "live.births"
}

cols <- c(y.var,"eid","agePulse","agePulse.squared","has.wes",paste0("PC",seq(1,40,by=1)))

phen.matrix <- UKBB.phenotype.data[sexPulse == sex & !is.na(get(y.var)) & eid %in% samples.UKBB.200k[,eid],..cols]
phen.matrix[,has.children:=if_else(get(y.var)>0,1,0)]

# Now restrict to individuals we want to test:
gt.dt <- gt.dt[eid %in% phen.matrix[,eid]]

# Make fam file:
gt.fam <- phen.matrix[,c("eid")]
setnames(gt.fam,"eid","FID")
gt.fam[,IID:=FID]
gt.fam[,father:=0]
gt.fam[,mother:=0]
gt.fam[,sex.code:=0]
gt.fam[,covar:=-9]
setkey(gt.fam,"FID")
# fwrite(gt.fam,paste0(outfile,".fam"),quote = F, col.names = F,sep=" ")

# Make MAP file:
genes <- genes.hg38[,c("GENE","CHR","START")]
genes[,cM:=0]
setcolorder(genes,c("CHR","GENE","cM","START"))
# fwrite(genes,paste0(outfile,".map"),quote = F, col.names = F, sep = "\t")

# Write sparse GT file (for editing later):
fwrite(gt.dt,paste0(outfile,".sparse.ped"),quote = F, col.names = F, sep = "\t")

# Make Covar file
cov.matrix <- phen.matrix[,c("eid","agePulse","agePulse.squared","has.wes",paste0("PC",seq(1,40)))]
setkey(cov.matrix,"eid")
setnames(cov.matrix,"eid","FID")
cov.matrix[,IID:=FID]
setcolorder(cov.matrix,c("FID","IID","agePulse","agePulse.squared","has.wes",paste0("PC",seq(1,40))))
# fwrite(cov.matrix,paste0(outfile,".covariates.txt"), sep = " ", quote = F)

# Make Pheno file
phenotype.matrix <- phen.matrix[,c("eid","has.children")]
setkey(phenotype.matrix,"eid")
setnames(phenotype.matrix,"eid","FID")
phenotype.matrix[,IID:=FID]
setcolorder(phenotype.matrix,c("FID","IID","has.children"))
# fwrite(phenotype.matrix,paste0(outfile,".phenotypes.txt"), sep = " ", quote = F)