#' @title vuong_test
#' @description
#' Vuong's test
#'
#' This function performs Vuong's test, a likelihood ratio test for model selection and
#' non-nested hypotheses.
#' This function is for model selection between zero-inflated Poisson model and zero-inflated negative binomial model.
#'
#' @param phenodata a data frame containing family and individual IDs for all objects as well as
#' zero-inflated counts as a phenotype and a set of covariates.
#' Each row represents a different individual.
#' The first two columns are Family ID (FID) and Individual ID (IID) respectively.
#' There must be one and only one phenotype in the third column and
#' the phenotype have to be zero-inflated count data which should be non-negative integers, e.g. neuritic plaque counts.
#' Each of the rest of columns represents a different covariate, e.g. age, sex, etc.
#' @returns nothing returned, prints a table of 3 test statistics and p values, and exits silently.
#' @importFrom pscl vuong
#' @export

vuong_test <- function(phenodata){
  dat <- phenodata[,-(1:2)]
  colnames(dat)[1] <- "y"
  zip <- zeroinfl(y~.|.,data = dat,EM=T)
  zinb <- zeroinfl(y~.|.,data = dat,dist = "negbin",EM=T)
  vuong(zip,zinb)
}

