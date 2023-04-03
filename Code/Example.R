# Set up enviornment
library(MGLM)
library(SKAT)
library(qvalue)
library(fda)
library(MASS)
library(survival)
library(globaltest)

source("GATI_MFG_function.R")

setwd("/Users/huanshan/Documents/Projects/Gene_by_gene/GenEpi_revision/GitHub")
# Phenotype: CHD status of infants
baby_pheno <- read.table("Fetal_phenotype.txt", header = T)
baby_pheno <- as.matrix(baby_pheno[,-1])

# Maternal genotype (100 mothers; 50 SNPs in the gene region)
mom_geno <- read.table("Maternal_genotype.txt", header = T)
mom_geno <- as.matrix(mom_geno [,-1])

# Maternal genotype (100 infants that paired with the 100 mothers; 50 SNPs in the gene region)
baby_geno <- read.table("Fetal_genotype.txt", header = T)
baby_geno <- as.matrix(baby_geno [,-1])

# Covariates
covariates <- read.table("Covariates.txt", header = T)
covariates <- as.matrix(covariates[,-1])

# Weights for SNPs in maternal gene, weights for SNPs in fetal gene, weights for the interaction 
weights_mom <- rep(1, ncol(mom_geno))
weights_baby <- rep(1, ncol(baby_geno))
weights_inte <- rep(1, ncol(mom_geno)*ncol(baby_geno))

# Evaluate the joint effect of maternal gene, fetal gene and their interaction using GATI-MFG
GATI_MFG(baby_pheno, mom_geno, baby_geno, covariates, weights_mom, weights_baby, weights_inte)



