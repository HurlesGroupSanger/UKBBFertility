#!/software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript --vanilla

library(GMMAT)
library(data.table)
library(tidyr)
library(dplyr)
library(Matrix)

args <- commandArgs(trailingOnly = T)

variant.counts <- readRDS("/lustre/scratch115/teams/hurles/users/eg15/UKBB/per_gene_tests/rawdata/snvresources/variant_counts.rdat")
phenotypes <- readRDS("/lustre/scratch115/teams/hurles/users/eg15/UKBB/per_gene_tests/rawdata/phenofiles/UKBB.phenotype.rdat")
grm.matrix <- readRDS("/lustre/scratch115/teams/hurles/users/eg15/UKBB/ukbb_penetrance/IBD/grm.MARTIN_IDS.sparsematrix.rds")

grm.allowed <- row.names(grm.matrix)

lm.table <- merge(phenotypes,variant.counts[type == "LOF_HC" & allele.freq == 0,c("sample_id","product_sHET")],by.x="eid",by.y="sample_id")
lm.table <- lm.table[eid %in% grm.allowed]
lm.table <- lm.table[sexPulse == 1]
lm.table <- lm.table[!is.na(children.fathered)]
lm.table[,has.children:=if_else(children.fathered>0,1,0)]

lm.table <- sample_n(lm.table,20)

samples <- lm.table[,eid]
lm.matrix <- grm.matrix[samples,samples]

covariates <- c("agePulse.squared","agePulse")
cov.string <- paste(covariates, collapse=" + ")
formated.formula <- as.formula(paste("has.children", cov.string,sep=" ~ "))

model0 <- glmmkin(formated.formula, data = lm.table, kins = lm.matrix, id = "eid", family = binomial(link = "logit"))
