#!/software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript --vanilla

library(data.table)
library(tidyr)
library(dplyr)
library(broom)
library(meta)

args <- commandArgs(trailingOnly = T)

input.icd <- args[1]
variant.maf <- as.integer(args[2])
age.cutoff <- as.integer(args[3])
fi.only <- as.logical(args[4])
outdir <- args[5]

disease.data.long <- readRDS(input.icd)
disease.data.long <- disease.data.long[age.at.incidence < age.cutoff]

icd.codes <- fread("rawdata/phewas/icd10_tree.tsv", header = T)
variant.counts <- readRDS("rawdata/phewas/variant_counts.rdat")
UKBB.phenotype.data <- readRDS("rawdata/phewas/UKBB.phenotype.rdat")
if (fi.only == T) {
  UKBB.phenotype.data <- UKBB.phenotype.data[has.first.incidence.data == T]
}

current.code <- as.integer(Sys.getenv("LSB_JOBINDEX"))
current.code.line <- icd.codes[current.code]

outfile <- paste("glm",current.code,"out",sep=".")
outfile <- paste(outdir,outfile,sep="/")

get.participants <- function(code, node, ids = c()) {
  
  ids <- c(ids, unique(disease.data.long[icd.code %in% code,eid]))
  curr.children <- icd.codes[parent_id == node]
  
  if (nrow(curr.children) != 0) {
    for (i in 1:nrow(curr.children)) {
      ids <- c(get.participants(curr.children[i, coding], curr.children[i, node_id], ids))
    }
    return(ids)
  } else {
    return(ids)
  }
}

run.lm <- function(code, node, sex, variant.type) {
  
  ## Some individuals have multiple codes in the same category when dealing with blocks
  indv.with.code <- unique(get.participants(code, node))
  
  lm.table <- merge(UKBB.phenotype.data,variant.counts[allele.freq == variant.maf & type == variant.type],by.x="eid",by.y="sample_id")
  
  lm.table[,has.disorder:=if_else(eid %in% indv.with.code,1,0)]
  
  if (variant.type == "DEL") {
    lm.table[has.wes == F]
  }
  
  if (sex == "MALE") {
    lm.table <- lm.table[sexPulse == 1 & !is.na(children.fathered)]
    lm.table[,childless:=if_else(children.fathered > 0,1,0)]
  } else if (sex == "FEMALE") {
    lm.table <- lm.table[sexPulse == 2 & !is.na(live.births)]
    lm.table[,childless:=if_else(live.births > 0,1,0)]
  }
  
  if (nrow(lm.table[has.disorder == 1]) < 2) {
    return(list(NaN,
                NaN,
                NaN,
                NaN,
                NaN,
                NaN,
                nrow(lm.table[has.disorder == 1]),
                nrow(lm.table[has.disorder == 1 & product_sHET > 0])))
  } else {
  
    covariates <- c("product_sHET","has.disorder","sexPulse","agePulse","agePulse.squared","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30")
    lm.formula <- as.formula(paste("childless",paste(covariates,collapse = "+"),sep="~"))
    
    lm.model <- glm(lm.formula, data = lm.table, family = "binomial")
    lm.out <- data.table(tidy(lm.model))

    return(list(lm.out[term=="product_sHET",estimate],
                lm.out[term=="product_sHET",std.error],
                lm.out[term=="product_sHET",p.value],
                lm.out[term=="has.disorder",estimate],
                lm.out[term=="has.disorder",std.error],
                lm.out[term=="has.disorder",p.value],
                nrow(lm.table[has.disorder == 1]),
                nrow(lm.table[has.disorder == 1 & product_sHET > 0])))
    
  }
  
}

analysis.table <- data.table(crossing(current.code.line[,c("coding","meaning","node_id")],sex=c("MALE","FEMALE"),variant.type = c("DEL","LOF_HC")))
analysis.table[,c("var.est","var.err","var.p","icd.est","icd.err","icd.p","n","n.cases"):=run.lm(coding,node_id,sex,variant.type),by=1:nrow(analysis.table)]

get.meta.val <- function(s) {
  
  meta.table <- analysis.table[sex == s]
  
  result <- list()
  
  if (nrow(meta.table[is.na(var.est)]) > 0) {
    
    result <- append(result, list(NaN,NaN,NaN))
    
  } else {
    
    meta.analy.var <- metagen(var.est,
                              var.err,
                              studlab = variant.type,
                              method.tau = "SJ",
                              sm = "OR",
                              data = meta.table)
    
    result <- append(result, list(meta.analy.var$TE.fixed,
                                  meta.analy.var$seTE.fixed,
                                  meta.analy.var$pval.fixed))
    
  }
  
  if (nrow(meta.table[is.na(icd.est)]) > 0) {
    
    result <- append(result, list(NaN,NaN,NaN))
    
  } else {
    
    meta.analy.icd <- metagen(icd.est,
                              icd.err,
                              studlab = variant.type,
                              method.tau = "SJ",
                              sm = "OR",
                              data = meta.table)
    
    result <- append(result, list(meta.analy.icd$TE.fixed,
                                  meta.analy.icd$seTE.fixed,
                                  meta.analy.icd$pval.fixed))
    
  }
    
  result <- append(result, list(meta.table[,sum(n)],meta.table[,sum(n.cases)]))
  
  return(result)
  
}

meta.table <- data.table(current.code.line[,c("coding","meaning","node_id")],
                         sex = c("MALE","FEMALE"),
                         variant.type = "META")

meta.table[,c("var.est","var.err","var.p","icd.est","icd.err","icd.p","n","n.cases"):=get.meta.val(sex),by=1:nrow(meta.table)]
analysis.table <- rbind(analysis.table, meta.table)

get.chapter <- function(code) {
  
  curr.code <- icd.codes[node_id == code,coding]
  curr.par <- icd.codes[node_id == code,parent_id]
  
  if (grepl("Chapter",curr.code) == T) {
    return(curr.code)
  } else if (grepl("-1",curr.code) == T) {
    return(NA)
  } else {
    get.chapter(curr.par)
  }
  
}

analysis.table[,chapter:=get.chapter(node_id),by=1:nrow(analysis.table)]

tree.climb <- function(code, level = 1) {
  
  curr.par <- icd.codes[node_id == code,parent_id]
  
  if (curr.par == 0) {
    return(level)
  } else {
    level <- tree.climb(curr.par, level + 1)
  }
  
}

analysis.table[,level:=tree.climb(node_id),by=1:nrow(analysis.table)]

write.table(analysis.table,outfile,col.names = F,row.names = F,sep = "\t", quote = F)






