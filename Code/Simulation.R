# Set up environment
library(MGLM)
library(SKAT)
library(qvalue)
library(fda)
library(MASS)
library(survival)
library(globaltest)

source("./FDA_functions/GFLM_fixed_model.R")
source("./FDA_functions/GFLM_beta_smooth_only.R")
source("GATI_MFG_function.R")
load("Haplotype_pool.RData")

ind_n <- 1000
snp_n <- 100
causal_snp_n <- 10

# simulation maternal and fetal genotypes
snp_index <- sample(ncol(hap_df), 300, replace = F)
ind_hap1 <- hap_df[sample(nrow(hap_df), ind_n, replace = T), snp_index]
ind_hap2 <- hap_df[sample(nrow(hap_df), ind_n, replace = T), snp_index]
ind_hap3 <- hap_df[sample(nrow(hap_df), ind_n, replace = T), snp_index]
genotype1 <- ind_hap1 + ind_hap2
genotype2 <- ind_hap1 + ind_hap3
genotype1[, colSums(genotype1)/(ind_n * 2) > 0.5] <- 2 - genotype1[, colSums(genotype1)/(ind_n * 2) > 0.5]
genotype2[, colSums(genotype2)/(ind_n * 2) > 0.5] <- 2 - genotype2[, colSums(genotype2)/(ind_n * 2) > 0.5]

polycol_index1 <- which(colSums(genotype1)  != 0)
polycol_index2 <- which(colSums(genotype2)  != 0)
polycol_index <- intersect(polycol_index1, polycol_index2)
loc <- sample(polycol_index, snp_n, replace = F)

genotype1 <- genotype1[, loc]
genotype2 <- genotype2[, loc]

gntp_int <- kr(genotype1, genotype2)

GT <- gntp_int/2
GT <- GT[, colSums(GT)  != 0] #remove nonpolymorphic

# specify multiplicative effect model
alpha <- 1.0
theta11 <- 0.2
theta12 <- 0.45
theta21 <- 0.2
theta22 <- 0.45

distmdl <- c(alpha, alpha * ( 1 + theta21), alpha * (1 + theta22),
             alpha * (1 + theta11), alpha * (1 + theta11) * (1 + theta21), alpha * (1 + theta11) * (1 + theta22),
             alpha * (1 + theta12), alpha * (1 + theta12) * (1 + theta21), alpha * (1 + theta12) * (1 + theta22))

names(distmdl) <- c('00', '01', '02', '10', '11', '12', '20', '21', '22')

# simulate phenotype
pheno <- runif(nrow(GT), 0, 1)

sn1 <- sample(ncol(genotype1), size = causal_snp_n, replace = F)
sn2 <- sample(ncol(genotype2), size = causal_snp_n, replace = F)

cb <- paste(genotype1[, sn1], genotype2[, sn2], sep = '')
vl <- matrix(distmdl[cb], ncol = causal_snp_n, nrow = ind_n)
eta <- rowSums(vl)
eta <- 0.52 * eta

eta <- eta - mean(eta) + log(0.47)
risk <- exp(eta)/(1 + exp(eta))
y <- as.numeric(pheno < risk)


# simulate covariates
x1 <- sample(c(0, 1), ind_n, replace = T, prob = c(0.5, 0.5))
x2 <- rnorm(ind_n)
covariates <- cbind(x1, x2)


#Evaluate overall effect of maternal gene and fetal gene using GATI-MFG
obj <- SKAT_Null_Model(y ~ covariates, out_type = "D")
psk1 <- SKAT(genotype1, obj, weights = rep(1, ncol(genotype1)))$p.value
psk2 <- SKAT(genotype2, obj, weights = rep(1, ncol(genotype2)))$p.value
psk3 <- SKAT(GT, obj, weights = rep(1, ncol(GT)))$p.value 

p_value_skat<- cauchy_combine(psk1, psk2, psk3)

# Evaluate overall effect of maternal gene and fetal gene using single-variant-based logistic regression (SLR)
comb <- expand.grid(1:ncol(genotype1), 1:ncol(genotype2)) #can use combn if # of snps same in genotype 1 and 2
pl <- rep(0, dim(comb)[1]) 
for(j in 1:length(pl)){
  combj <- as.numeric(comb[j, ])
  slm2 <- glm(y ~ x1 + x2 + genotype1[, combj[1]] + genotype2[, combj[2]] + genotype1[, combj[1]] * genotype2[, combj[2]], family = 'binomial')
  if (all(dim(summary(slm2)$coefficients) == c(6,4))){
    slm1 <- glm(y ~ x1 + x2, family = "binomial")
    mdlvs <- anova(slm1, slm2, test = "LRT")
    pl[j] <- mdlvs[2,5]
  }else{
    pl[j] <- NA
  }
}

ql <- qvalue(pl)$qvalues
rj <- min(ql, na.rm = T) < 0.05


## Evaluate overall effect of maternal gene and fetal gene using functional data analysis
pos1  <- c(1:ncol(genotype1))
pos2  <- c(1:ncol(genotype2))
pos3 <- c(1:ncol(GT))

order  =  4
bbasis = 10
gbasis = 10
fbasis = 11
gfasis = 11

gflm_fixed_bspline_pvalue1  = gflm_fixed_model(y, mode = "Additive", genotype1, pos1, order, bbasis, fbasis, gbasis, as.matrix(covariates), base = "bspline", interaction = FALSE)
fda_p1 <- gflm_fixed_bspline_pvalue1$Rao

gflm_fixed_bspline_pvalue1  = gflm_fixed_model(y, mode = "Additive", genotype2, pos2, order, bbasis, fbasis, gbasis, as.matrix(covariates), base = "bspline", interaction = FALSE)
fda_p2 <- gflm_fixed_bspline_pvalue1$Rao

fml_test <- gflm_fixed_model(y, mode = "Additive", GT, pos3, order, bbasis, fbasis, gbasis, as.matrix(covariates), base = "bspline", interaction = FALSE)
fml <- fml_test$Rao


Pval_fda <- c(fda_p1, fda_p2, fml)
stat_fda <- mean(tan((0.5-Pval_fda)*pi))
p_value_fda <- 0.5-atan(stat_fda)/pi



# Put the p-values by three methods together and save the results
p_results <- c(p_value_t1e, rje_0.05, rje_0.01, rje_0.001, p_value_t1e_fda)
write.table(p_results, file = "simulation_multiplicative_results")














