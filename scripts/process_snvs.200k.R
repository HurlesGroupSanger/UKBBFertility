#!/software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript --vanilla

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

name.mod <- ""

## Load gene lists:
## Load gene translation file
gene.translate <- fread("hg38_to_hg19_ENSG.txt", header=F)
setnames(gene.translate,names(gene.translate),c("hg38.GENE","hg19.GENE"))

## Load sHET genes
shet.genes <- fread("shet.hgnc.txt")
setnames(shet.genes,names(shet.genes),c("hg19.GENE","GENE","sHET.val","HGNC.ID"))
shet.genes[,deciles:=cut(sHET.val,breaks=quantile(sHET.val,seq(0,1,by=0.1)),include.lowest = T)]
shet.genes[,sHET.val.binary:=cut(sHET.val,breaks=c(0,0.15,1),labels = c("lt_015","gt_015"),right = F)]

## Load Cassa sHET genes
shet.genes.cassa <- fread("shet.cassa.hgnc.txt")
## Remember, we had to annotate sHET with hg38 gene IDs
setnames(shet.genes.cassa,names(shet.genes.cassa),c("hg19.GENE","GENE","sHET.val","HGNC.ID"))

## Load pLI genes:
pli.genes <- fread("hg19.all_genes_with_pli.txt")
pli.genes <- pli.genes[,c("V1","V6","V8","V5")]
setnames(pli.genes, c("hg19.GENE","GENE","pLI.val","HGNC.ID"))
pli.genes[,pLI.val.binary:=cut(pLI.val,breaks=c(0,0.9,1),labels = c("lt_09","gt_09"),right = F)]

## Generate lists for phenotype testing:
gene.lists <- list()
gene.lists["highPLI"] = list(pli.genes[pLI.val >= 0.9,hg19.GENE])
gene.lists["highsHET"] = list(shet.genes[sHET.val >= 0.15,hg19.GENE])

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
  
  ## Filter by PEXT:
  # Only filter by PEXT for allosomes -- hemizygous genes do not use PEXT
  # variant.counts <- variant.counts[PEXT >= 0.1]
  
  ## Filter by allele missingness:
  variant.counts <- variant.counts[AN >= AP/2]
  ###!! FILTERING !!###
  
  return(list(variant.counts,final.samples))
  
}

# UKBB/Hg38
counts.file <- "counts.ukbb_wes.200k.txt"
samples.file <- "ukbb_ids.200k.txt"
indv.to.keep <- fread("indv_to_keep.txt",header = F)
indv.to.keep <- indv.to.keep[,V1]

## The weird paste is because of how sampleIDs are represented in the variants table
data.UKBB.200k <- load.data(counts.file,samples.file,indv.to.keep)
variants.UKBB.200k <- data.UKBB.200k[[1]]
samples.UKBB.200k <- data.UKBB.200k[[2]]
rm(data.UKBB.200k)

## Translate geneIDs from Hg38 to Hg19
variants.UKBB.200k <- merge(variants.UKBB.200k,gene.translate,by.x = "GENE",by.y = "hg38.GENE",all.x=T)
variants.UKBB.200k[,GENE:=hg19.GENE]
variants.UKBB.200k[,hg19.GENE:=NULL]

## Count Variants
count.variants <- function(variant.counts,samples,af.cutoff) {
  
  counts.AC1 <- variant.counts[UKBB.AC == 1,sum(gt),by=c("eid","CSQ")]
  counts.AC1[,AF:="AC1"]
  counts.AF0.1 <- variant.counts[UKBB.AF < af.cutoff,sum(gt),by=c("eid","CSQ")]
  counts.AF0.1[,AF:="AF0.1"]
  
  count.table <- bind_rows(counts.AC1,counts.AF0.1)
  
  samples.table <- data.table(crossing(samples=samples[,eid],CSQ=c("LOF_HC","MIS","SYN"),AF=c("AC1","AF0.1")))
  setnames(samples.table,"samples","eid")
  
  count.table <- merge(samples.table,count.table,by=c("eid","CSQ","AF"),all.x=T)
  count.table[,V1:=ifelse(is.na(V1),0,V1)]
  
  setnames(count.table,"V1","count")
  
  return(count.table)
  
}

UKBB.counts <- count.variants(variants.UKBB.200k,samples.UKBB.200k,0.001)
saveRDS(UKBB.counts,paste0(name.mod,"counts.rds"))

calc.gene.counts <- function(freq) {

  if (freq == 0) {
    UKBB.genes.200k <- unique(variants.UKBB.200k[UKBB.AC == 1,c("chr","pos","ref","alt","CSQ","GENE","gt")])
  } else {
    UKBB.genes.200k <- unique(variants.UKBB.200k[UKBB.AF < freq,c("chr","pos","ref","alt","CSQ","GENE","gt")])
  }
  UKBB.genes.200k <- UKBB.genes.200k[,sum(gt),by=c("GENE","CSQ")]
  setnames(UKBB.genes.200k,"V1","UKBB")
  UKBB.genes.200k[,maf:=freq]
  return(UKBB.genes.200k)
    
}

UKBB.genes.200k <- bind_rows(calc.gene.counts(0.001),
                             calc.gene.counts(0.0001),
                             calc.gene.counts(0.00001),
                             calc.gene.counts(0))
saveRDS(UKBB.genes.200k,paste0(name.mod,"genes.rds"))

rm(UKBB.genes.200k, UKBB.counts)

## Get per gene counts:
get.gene.counts.snvs <- function(data, samples, maf) {
  
  variants <- data.table(crossing(sample_id=samples[,eid],type=c("LOF_HC","SYN","MIS")))
  
  ## This handles the high pLI/sHET lists:
  for (gene.list in names(gene.lists)) {
    
    current.genes <- get(gene.list,gene.lists)
    
    ## Filter on frequency cutoff and gene list
    if (maf == 0) {
      counts <- data[UKBB.AC == 1 & GENE %in% current.genes & chr!="X" & chr!="Y",sum(gt),by=c("eid","CSQ")]
    } else {
      counts <- data[UKBB.AF <= maf & GENE %in% current.genes & chr!="X" & chr!="Y",sum(gt),by=c("eid","CSQ")]
    }
    
    variants <- merge(variants,counts,by.x=c("sample_id","type"),by.y=c("eid","CSQ"),all.x=T)
    variants[,V1:=if_else(is.na(V1),0L,V1)]
    setnames(variants,"V1",gene.list)
    
  }
  
  variants[,allele.freq:=maf]
  setcolorder(variants,c(names(variants)[-grep("type",names(variants))],c("type")))
  return(variants)
}

snv.counts <- bind_rows(get.gene.counts.snvs(variants.UKBB.200k,samples.UKBB.200k,0),
                        get.gene.counts.snvs(variants.UKBB.200k,samples.UKBB.200k,1e-5),
                        get.gene.counts.snvs(variants.UKBB.200k,samples.UKBB.200k,1e-4),
                        get.gene.counts.snvs(variants.UKBB.200k,samples.UKBB.200k,1e-3))

## Now do quantitative calculations:
quantify.shet.snv <- function(maf, type, list.name, shet.to.use, to.remove=c()) {
  
  ## Filter on frequency cutoff and gene list
  if (maf == 0) {
    returned <- variants.UKBB.200k[UKBB.AC == 1 & chr!="X" & chr!="Y" & CSQ == type & !GENE %in% to.remove]
  } else {
    returned <- variants.UKBB.200k[UKBB.AF <= maf & chr!="X" & chr!="Y" & CSQ == type & !GENE %in% to.remove]
  }
  
  returned <- merge(returned,shet.to.use[,c("hg19.GENE","sHET.val")],by.x="GENE",by.y="hg19.GENE")
  returned <- returned[,list(1-prod(1-sHET.val)),by="eid"]
  
  setnames(returned,c("eid","V1"),c("sample_id",list.name))
  returned[,allele.freq:=maf]
  returned[,type:=type]
  
  return(returned)
  
}

## load disease genes
disease.genes <- fread("diseaseGenes.txt",header=F)
setnames(disease.genes,"V1","hg19.GENE")

## Load male infertility genes
male.infertility.genes <- fread("male_infertility_genes.annotated.txt",header=F)
setnames(male.infertility.genes,names(male.infertility.genes),c("GENE","inheritance","evidence","hg19.GENE"))

## Load mouse infertility genes
mouse.infertility.genes <- fread("mouse_infertility_genes.annotated.txt",header=F)
setnames(mouse.infertility.genes,names(mouse.infertility.genes),c("GENE","hg19.GENE","chr","start"))

## Make lists to make the below script work faster
exclusion.list <- list("product_sHET"=c(),
                       "product_sHET_old"=c(),
                       "product_sHET_no_maleInfertilityGenes" = male.infertility.genes[,hg19.GENE],
                       "product_sHET_no_diseaseGenes" = disease.genes[,hg19.GENE],
                       "product_sHET_no_mouseInfertilityGenes" = mouse.infertility.genes[,hg19.GENE])

## Various expression cutoffs:
files <- list.files("tissues/")

for (tissue in files) {
  
  t.dt <- fread(paste0("tissues/",tissue), header = F)
  setnames(t.dt,"V1","hg19.GENE")
  t.name <- str_remove(tissue,".txt")
  t.list <- paste0("product_sHET_no_", t.name)
  exclusion.list[[t.list]] <- t.dt[,hg19.GENE]
  
}

shet.quant.table <- data.table(crossing(allele.freq=unique(snv.counts[,allele.freq]),type=unique(snv.counts[,type])))

for (set in names(exclusion.list)) {
  
  quants <- data.table()
  
  for (i in c(1:nrow(shet.quant.table))) {
    
    if (set == "product_sHET_old"){
      shet.to.use <- shet.genes.cassa    
    } else {
      shet.to.use <- shet.genes
    }
    
    returned <- quantify.shet.snv(shet.quant.table[i,allele.freq],
                                  shet.quant.table[i,type],
                                  set,
                                  shet.to.use,
                                  to.remove = get(set,exclusion.list))
    
    returned[,allele.freq:=shet.quant.table[i,allele.freq]]
    returned[,type:=shet.quant.table[i,type]]
    
    quants <- bind_rows(quants,returned)
    
  }
  
  snv.counts <- merge(snv.counts,quants,by=c("sample_id","allele.freq","type"),all.x=T)
  snv.counts[,eval(set):=if_else(is.na(get(set)),0,get(set))]
  
}

saveRDS(snv.counts,paste0(name.mod,"200k_counts.rds"))

rm(quants,shet.quant.table,returned)