#!/software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript --vanilla

library(data.table)
library(tidyr)
library(dplyr)
library(broom)
library(rcompanion)

## Load Necessary Data
UKBB.phenotype.data <- readRDS("rawdata/UKBB.phenotype.rdat")
variant.counts <-  readRDS("rawdata/variant_counts.rdat")
lm.master.table <- readRDS("rawdata/lm.master.table.rdat")

## Pick the right row to do from lm.master.table:
current.code <- as.integer(Sys.getenv("LSB_JOBINDEX"))

current.line <- lm.master.table[current.code]

run.regression <- function(maf,
                           gene.list,
                           y.var,
                           sex,
                           variant.type,
                           model.family,
                           add.covars = c(),
                           return.data = F,
                           remove.zeros = F,
                           remove.sequenced = T,
                           cutoff.high = F,
                           num.pcs = 40,
                           num.rare.pcs = 100,
                           indv.to.exclude = "",
                           indv.to.exclude.value = -1,
                           nagel = F) {
  
  id.name <- "eid"
  
  cols.to.keep <- c("sample_id",gene.list)
  final.stats <- variant.counts[type == variant.type & allele.freq == maf,..cols.to.keep]
  
  if (cutoff.high == T & (variant.type == "DEL" | variant.type == "DUP" | variant.type == "LOF_HC")) {
    final.stats <- final.stats[get(gene.list) <= 3]
  }
  
  final.stats <- merge(final.stats,UKBB.phenotype.data,by.x="sample_id",by.y=id.name)
  
  ## Remove WES individuals from CNV analyses for meta-analysis purposes
  if (remove.sequenced == T & (variant.type == "DEL" | variant.type == "DUP")) {
    final.stats <- final.stats[has.wes == 0]
  } else {
    add.covars <- c(add.covars,"has.wes")
  }
  
  ## Remove missing y.var data
  final.stats <- final.stats[!is.na(get(y.var))]
  if (remove.zeros == T) {
    final.stats <- final.stats[get(y.var)>0]
  }
  
  ## Do a model w/o sex
  if (sex == 1 | sex == 2) {
    final.stats <- final.stats[sexPulse == sex]
  }
  
  ## Remove individuals based on some criteria passed to this function
  if (indv.to.exclude != "") {
    to.exclude.list <- UKBB.phenotype.data[get(indv.to.exclude) == indv.to.exclude.value, eid]
    final.stats <- final.stats[!sample_id %in% to.exclude.list]  
  }
  
  ## Set linear or logistic model
  if (model.family == "binomial") {
    ## And force the phenotype to binary:
    final.stats[,binary.stat:=if_else(get(y.var) > 0,1,0)]
    y.var <- "binary.stat"
  }
  
  ## Remove missing additional covar data
  add.covars <- unlist(add.covars)
  for (cov in add.covars) {
    final.stats <- final.stats[!is.na(get(cov))]
  }
  
  covariates <- c(gene.list,add.covars,"agePulse.squared","agePulse","sexPulse","birth.year.cut")
  ## Set numbers of common/rare ancestry PCs
  if (num.pcs == 40) {
    covariates <- c(covariates,paste0("PC",seq(1,40)))
  } else if (num.pcs > 0) {
    for (pc.num in c(1:num.pcs)) {
      covariates <- c(covariates, paste0("PC",pc.num))
    }
  }
  if (num.rare.pcs == 100) {
    covariates <- c(covariates,paste0("scaled.rare.PC",seq(1,100)))
  } else if (num.rare.pcs > 0) {
    for (pc.num in c(1:num.rare.pcs)) {
      covariates <- c(covariates, paste0("scaled.rare.PC",pc.num))
    }
  }
  
  cov.string <- paste(covariates, collapse=" + ")
  formated.formula <- as.formula(paste(y.var, cov.string,sep=" ~ "))
  
  test.lm <- glm(formated.formula, data=final.stats, family=model.family)
  coef.lm <- tidy(test.lm) %>% data.table()
  final.stats <- augment(test.lm) %>% data.table()
  total.hits <- final.stats[,sum(get(gene.list))]
  
  ## NAGLEKERKE STUFF
  if (nagel == T) {
    covariates <- c("agePulse.squared","agePulse","sexPulse","birth.year.cut",add.covars,paste0("PC",seq(1,40)),paste0("scaled.rare.PC",seq(1,100)))
    
    cov.string <- paste(covariates, collapse=" + ")
    formated.formula <- as.formula(paste(y.var, cov.string,sep=" ~ "))
    
    test.lm.init <- glm(formated.formula, data=final.stats, family=model.family)
    nag.shet.out <- nagelkerke(test.lm, null = test.lm.init)
    
    ## NAGLEKERKE STUFF
    if (return.data == T) {  
      return(list(coef.lm[term==eval(gene.list),estimate],
                  coef.lm[term==eval(gene.list),std.error],
                  coef.lm[term==eval(gene.list),p.value],
                  total.hits,
                  nrow(final.stats),
                  list(coef.lm),
                  nag.shet.out$Pseudo.R.squared.for.model.vs.null[3]))
    } else {
      return(list(coef.lm[term==eval(gene.list),estimate],
                  coef.lm[term==eval(gene.list),std.error],
                  coef.lm[term==eval(gene.list),p.value],
                  total.hits,
                  nrow(final.stats),
                  list(),
                  nag.shet.out$Pseudo.R.squared.for.model.vs.null[3]))
    }
    
  } else {
    
    if (return.data == T) {
      return(list(coef.lm[term==eval(gene.list),estimate],
                  coef.lm[term==eval(gene.list),std.error],
                  coef.lm[term==eval(gene.list),p.value],
                  total.hits,
                  nrow(final.stats),
                  list(coef.lm),
                  NaN))
    } else {
      return(list(coef.lm[term==eval(gene.list),estimate],
                  coef.lm[term==eval(gene.list),std.error],
                  coef.lm[term==eval(gene.list),p.value],
                  total.hits,
                  nrow(final.stats),
                  list(),
                  NaN))
    }
    
  }
  
}

result <- current.line[,c("var.beta","var.stderr","var.p","n.var","n.indvs","model","inc.r.shet"):=run.regression(maf,
                                                                                                                  gene.list,
                                                                                                                  y.var,
                                                                                                                  sex,
                                                                                                                  variant.type,
                                                                                                                  model.family,
                                                                                                                  add.covars = add.covars,
                                                                                                                  return.data = return.data,
                                                                                                                  remove.zeros = remove.zeros,
                                                                                                                  remove.sequenced = remove.sequenced,
                                                                                                                  cutoff.high = cutoff.high,
                                                                                                                  num.pcs = num.pcs,
                                                                                                                  num.rare.pcs = num.rare.pcs,
                                                                                                                  indv.to.exclude = indv.to.exclude,
                                                                                                                  indv.to.exclude.value = indv.to.exclude.value,
                                                                                                                  nagel = nagel)]

file <- paste0("outfiles/result.",current.code,".rdat")
saveRDS(result, file)



