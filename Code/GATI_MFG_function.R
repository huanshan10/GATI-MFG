GATI_MFG <- function(baby_pheno, mom_geno, baby_geno, covariates, weights_mom, weights_baby, weights_inte){
  # interaction genotype
  gntp_int <- kr(mom_geno, baby_geno)
  GT <- gntp_int/2
  GT <- GT[, colSums(GT)  != 0] #remove nonpolymorphic
  weights_inte2 <- weights_inte[colSums(GT)  != 0]
  
  # test for the main effect of maternal gene, main effect of fetal gene and interaction effect individually
  obj_null <- SKAT_Null_Model(baby_pheno ~ covariates, out_type = "D")
  psk1 <- SKAT(mom_geno, obj_null, weights = weights_mom)$p.value
  psk2 <- SKAT(baby_geno, obj_null, weights = weights_baby)$p.value
  psk3 <- SKAT(GT, obj_null, weights = weights_inte2)$p.value  
  
  # combine the three p-values using Cauchy combination test
  p_value_overall<- cauchy_combine(psk1, psk2, psk3)
  
  return(p_value_overall)
}


#function to combine p-values using Cauchy combination test
cauchy_combine <- function(x1, x2, x3){
  x <- c(x1, x2, x3)
  stat <- mean(tan((0.5-x)*pi))
  p_value <- 0.5-atan(stat)/pi
  return(p_value)
}
