#!/software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript --vanilla

library(data.table)
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly = T)

variant.maf <- as.integer(args[1])
variant.type <- args[2]
sex <- as.integer(args[3])
outfile <- args[4]

# Load the phenotype matrix:
UKBB.phenotype.data <- readRDS("rawdata/phenofiles/UKBB.phenotype.rdat")
# Load the variant counts:
variant.counts <- readRDS("rawdata/snvresources/variant_counts.rdat")

# Building covariate matrix:
if (sex == 1) {
  y.var <- "children.fathered"
} else {
  y.var <- "live.births"
}

cols <- c(y.var,"eid","agePulse","agePulse.squared","has.wes",paste0("PC",seq(1,40,by=1)))

phen.matrix <- UKBB.phenotype.data[sexPulse == sex & !is.na(get(y.var)) & has.wes > 0,..cols]
phen.matrix[,has.children:=if_else(get(y.var)>0,1,0)]

# Now restrict to individuals we want to test:
phen.matrix <- merge(phen.matrix, variant.counts[type == variant.type & allele.freq == variant.maf & sample_id %in% phen.matrix[,eid],c("sample_id","product_sHET")],by.x = "eid",by.y = "sample_id")
setkey(phen.matrix,"eid")
phen.matrix[,highsHET:=if_else(product_sHET>0.15,"12","11")]

# Make fam file:
gt.fam <- phen.matrix[,c("eid","highsHET")]
setnames(gt.fam,"eid","FID")
gt.fam[,IID:=FID]
gt.fam[,father:=0]
gt.fam[,mother:=0]
gt.fam[,sex.code:=0]
gt.fam[,covar:=-9]
setkey(gt.fam,"FID")
fwrite(gt.fam[,c("FID","IID","father","mother","sex.code","covar")],paste0(outfile,".fam"),quote = F, col.names = F,sep=" ")

# Write PED file:
fwrite(gt.fam[,c("FID","IID","father","mother","sex.code","covar","highsHET")],paste0(outfile,".ped"),quote = F, col.names = F, sep = "\t")

# Make MAP file:
genes <- data.table("CHR" = 1, "GENE" = "SHET", "cM"= 0, "START" = 10)
fwrite(genes,paste0(outfile,".map"),quote = F, col.names = F, sep = "\t")

# Make Covar file
cov.matrix <- phen.matrix[,c("eid","agePulse","agePulse.squared","has.wes",paste0("PC",seq(1,40)))]
setkey(cov.matrix,"eid")
setnames(cov.matrix,"eid","FID")
cov.matrix[,IID:=FID]
setcolorder(cov.matrix,c("FID","IID","agePulse","agePulse.squared","has.wes",paste0("PC",seq(1,40))))
fwrite(cov.matrix,paste0(outfile,".covariates.txt"), sep = " ", quote = F)

# Make Pheno file
phenotype.matrix <- phen.matrix[,c("eid","has.children")]
setkey(phenotype.matrix,"eid")
setnames(phenotype.matrix,"eid","FID")
phenotype.matrix[,IID:=FID]
setcolorder(phenotype.matrix,c("FID","IID","has.children"))
fwrite(phenotype.matrix,paste0(outfile,".phenotypes.txt"), sep = " ", quote = F)