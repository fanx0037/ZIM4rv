#' Estimation of phi_hat and lambda_hat for ZIP model
#'
#' This function gives the estimation of 2 parameters phi and lambda for each subject
#' under the null hypothesis.
#'
#' @param simud a data frame containing a phenotype named y and covariates
#' @return a list of 2 estimations of parameters for each subject
#' @details
#' This function first fits zero‐inflated Poisson regression of phenotype y
#' on the covariates only to obtain the estimates of regression coefficients
#' and then compute the estimations of phi and lambda.
#' @seealso \code{zeroinfl}
#' @importFrom pscl zeroinfl
#' @importFrom CompQuadForm davies
#' @importFrom stats cov
#' @importFrom stats dbeta
#' @importFrom stats pchisq
#' @export


phi_lambda_hat <- function(simud){
  simudata <- as.data.frame(simud)
  m1 <- zeroinfl(y ~ ., data = simudata)
  a1_hat <- m1$coefficients$zero
  a2_hat <- m1$coefficients$count
  X <- as.matrix(cbind(rep(1, nrow(simudata)), simudata[, -1]))
  hat_fi <- exp(X%*%a1_hat)
  hat_lambda <- exp(X%*%a2_hat)
  hat <- list(hat_fi = hat_fi, hat_lambda = hat_lambda)
  return(hat)
}


#' Compute Score statistics for ZIP model
#'
#' This function takes the estimations of phi and lambda produced by the \code{phi_lambda_hat}
#' and computes the score statistics under the null hypothesis.
#'
#' @param simudata a data frame containing a phenotype named y and covariates
#' @param G_rare a data frame containing data of rare variants with the same subject order as in simudata
#' @return a list of 2 matrice of the score statistics for each variant from each subject
#' @importFrom pscl zeroinfl
#' @importFrom CompQuadForm davies
#' @importFrom stats cov
#' @importFrom stats dbeta
#' @importFrom stats pchisq
#' @export

U_fi_lmd <- function(simudata, G_rare){
  hat_fi <- phi_lambda_hat(simudata)$hat_fi
  hat_lambda <- phi_lambda_hat(simudata)$hat_lambda
  comb <- cbind(simudata$y, hat_fi, hat_lambda)
  colnames(comb)[2:3] <- c( "fi","lambda")
  comb <- as.matrix(comb)
  u_fi <- matrix(nrow = nrow(G_rare), ncol = ncol(G_rare))

  for (i in 1:nrow(G_rare)){

    for (j in 1:ncol(G_rare)){

      if (comb[i,1]==0) {
        fim <- comb[i,2]*G_rare[i,j]/(comb[i,2] + exp(-comb[i,3])) - comb[i,2] * G_rare[i,j]/(1 + comb[i,2])
      }
      else {
        fim <- (- comb[i,2] * G_rare[i,j]/(1 + comb[i,2]))
      }
      u_fi[i,j] <- fim
    }
  }

  u_lambda=matrix(nrow=nrow(G_rare), ncol = ncol(G_rare))

  for (i in 1:nrow(G_rare)){

    for (j in 1:ncol(G_rare)){

      if (comb[i,1]==0) {
        fim= -exp(-comb[i,3]) * comb[i,3] * G_rare[i,j]/(comb[i,2] + exp(-comb[i,3])) }
      else { fim = (comb[i,1] - comb[i,3]) * G_rare[i,j]}

      u_lambda[i,j]=fim
    }
  }

  score=list(u_fi=u_fi, u_lambda=u_lambda)
  return(score)
}


#' Estimation of phi_hat, mu_hat and alpha_hat for ZINB model
#'
#' This function gives the estimation of 3 parameters phi, mu and alpha in ZINB model
#' for each subject under the null hypothesis.
#'
#' @param simud a data frame containing a phenotype named y and covariates
#' @return a list of 3 estimations of parameters for each subject
#' @details
#' This function first fits zero‐inflated negative binomial regression of phenotype y
#' on the covariates only to obtain the estimates of regression coefficients
#' and inverse dispersion
#' and then compute the estimations of phi, mu and alpha.
#' @seealso \code{zeroinfl}
#' @importFrom pscl zeroinfl
#' @importFrom CompQuadForm davies
#' @importFrom stats cov
#' @importFrom stats dbeta
#' @importFrom stats pchisq
#' @export


phi_mu_hat4zinb <- function(simud){
  simudata <- as.data.frame(simud)
  m1 <- zeroinfl(y ~ .|., data = simudata, dist = "negbin")
  a1_hat <- m1$coefficients$zero
  a2_hat <- m1$coefficients$count
  hat_alpha <- (m1$theta)^(-1)
  X <- as.matrix(cbind(rep(1, nrow(simudata)), simudata[, -1]))
  hat_phi <- exp(X%*%a1_hat)
  hat_mu <- exp(X%*%a2_hat)
  hat <- list(hat_phi = hat_phi, hat_mu = hat_mu, hat_alpha = hat_alpha)
  return(hat)
}



#' Compute score statistics for ZINB model
#'
#' This function takes the estimations of phi and lambda produced by the \code{phi_lambda_hat4negbin}
#' and computes the score statistics for ZINB model under the null hypothesis.
#'
#' @param simudata a data frame containing a phenotype named y and covariates
#' @param G_rare a data frame containing data of rare variants with the same subject order as in simudata
#' @return a list of 2 matrice of the score statistics for each variant from each subject
#' @importFrom pscl zeroinfl
#' @importFrom CompQuadForm davies
#' @importFrom stats cov
#' @importFrom stats dbeta
#' @importFrom stats pchisq
#' @export

U_phi_mu4zinb <- function(simudata, G_rare){
  estimation <- phi_mu_hat4zinb(simudata)
  hat_phi <- estimation$hat_phi
  hat_mu <- estimation$hat_mu
  hat_alpha <- estimation$hat_alpha
  comb <- cbind(simudata$y, hat_phi, hat_mu)
  colnames(comb)[2:3] <- c( "phi","mu")
  comb <- as.matrix(comb)

  u_phi <- matrix(nrow = nrow(G_rare), ncol = ncol(G_rare))

  for (i in 1:nrow(G_rare)){

    for (j in 1:ncol(G_rare)){

      if (comb[i,1]==0) {
        phi_im <- (1/(comb[i,2]+(1+hat_alpha*comb[i,3])^(-hat_alpha^(-1)))-1/(1+comb[i,2]))*comb[i,2]*G_rare[i,j]
      }
      else {
        phi_im <- -comb[i,2]*G_rare[i,j]/(1+comb[i,2])
      }
      u_phi[i,j] <- phi_im
    }
  }

  u_mu <- matrix(nrow=nrow(G_rare), ncol = ncol(G_rare))

  for (i in 1:nrow(G_rare)){

    for (j in 1:ncol(G_rare)){

      if (comb[i,1]==0) {
        mu_im <- -comb[i,3]*(1+hat_alpha*comb[i,3])^(-1-hat_alpha^(-1))*G_rare[i,j]/(comb[i,2]+(1+hat_alpha*comb[i,3])^(-hat_alpha^(-1)))
      }
      else {
        mu_im <- (comb[i,1]-comb[i,3])*G_rare[i,j]/(1+hat_alpha*comb[i,3])
      }

      u_mu[i,j] <- mu_im
    }
  }

  score=list(u_phi=u_phi, u_mu=u_mu)
  return(score)
}



#' Compute the p-value for the burden test
#'
#' This function takes a vector of weights, a data frame of rare variants
#' and a matrix of Score statistics produced by \code{U_fi_lmd} for ZIP model
#' or \code{U_phi_mu4zinb} for ZINB model to compute the p-value for the burden test.
#'
#' @param wt a numeric vector containing weights for all variants
#' @param G_rare a data frame containing data of rare variants
#' @param s a matrix of the score statistics for each variant from each subject
#' @return the p-value for the burden test
#' @importFrom pscl zeroinfl
#' @importFrom CompQuadForm davies
#' @importFrom stats cov
#' @importFrom stats dbeta
#' @importFrom stats pchisq
#' @export

p_burden_single = function(wt,G_rare, s){

  sum_var=rep(1,ncol(G_rare))
  for (i in 1:nrow(G_rare)){
    sum_var[i]=sum(wt*s[i,])
  }
  Q=(sum(sum_var)^2)/sum(sum_var^2)
  if (is.nan(Q)) {Q =0 }
  p=1-pchisq(Q,df=1)

  return(p)
}


#' Compute the p-value for the kernel test
#'
#' This function takes a diagonal matrix of weights, a data frame of rare variants
#' and a matrix of Score statistics produced by \code{U_fi_lmd} for ZIP model
#' or \code{U_phi_mu4zinb} for ZINB model to compute the p-value for the kernel test.
#'
#' @param wt_matrix2 a diagonal matrix containing the squared weights for all variants
#' @param G_rare a data frame containing data of rare variants
#' @param s a matrix of the score statistics for each variant from each subject
#' @return the p-value for the kernel test (ZIP-k)
#' @importFrom pscl zeroinfl
#' @importFrom CompQuadForm davies
#' @importFrom stats cov
#' @importFrom stats dbeta
#' @importFrom stats pchisq
#' @export

p_kernel_single = function(wt_matrix2,G_rare, s){

  onematrix=as.matrix(rep(1,nrow(G_rare)))
  T=as.vector(t(onematrix)%*%s%*% wt_matrix2 %*% t(s)%*%onematrix)
  V=nrow(G_rare)*cov(s)%*% wt_matrix2
  lambda_V=eigen(V, symmetric = T, only.values = T)$values
  p=davies(T, lambda = lambda_V)$Qq

  return(p)
}


#' Cauchy combination test (Cauchy‐p)
#'
#' This function combines p-values using Cauchy combination test
#' for the omnibus test of testing the joint genetic effect.
#'
#' @param x a numeric vector containing p-values
#' @return a combined p-value indicating the joint effect
#' @importFrom pscl zeroinfl
#' @importFrom CompQuadForm davies
#' @importFrom stats cov
#' @importFrom stats dbeta
#' @importFrom stats pchisq
#' @export

cauchyp<-function(x){
  pi=3.141593
  T<-0.5*tan((0.5-x[1])*pi) + 0.5*tan((0.5-x[2])*pi)
  p = 0.5 - (atan(T))/pi
  return(p)
}



#' Preprocess files in PLINK format
#'
#' This function converts PLINK format files into data frames containing genotypes,
#' phenotypes and covariates information in proper format.
#'
#' @param gene_name a character string of the name of a gene, e.g."CEPT". The name is case-sensitive.
#' @param region_file region file in PLINK format
#' @param dosage_file dosage file in PLINK format
#' @param fam_file .fam file in PLINK format
#' @param pheno_file phenotype file in PLINK format
#' @param cov_file covariate file in PLINK format
#' @returns a list of 2 data frames containing genotypes, phenotypes and covariates respectively
#' @export
#'
preprocess <- function(gene_name,region_file,dosage_file,fam_file,pheno_file,cov_file){
  # region + dosage + fam -> genedata

  # locate the region of the gene
  gene_ind <- which(region_file[,4]==gene_name)
  chr_ind <- region_file[gene_ind,1]
  gene_start <- region_file[gene_ind,2]
  gene_end <- region_file[gene_ind,3]
  # select the variants within the region of the gene
  dosage_sub <- dosage_file[dosage_file[,1]==chr_ind,]
  dosage_sub[dosage_sub=='na'] <- NA
  variants_ind <- which(dosage_sub[,4]>=gene_start & dosage_sub[,4]<=gene_end)
  dosage_transpose <- matrix(as.numeric(t(dosage_sub[variants_ind,-(1:6)])),ncol = length(variants_ind))
  genedata <- cbind(fam_file[,1:2],dosage_transpose)
  colnames(genedata) <- c("FID","IID")

  # pheno + cov -> phenodata
  pheno_file[pheno_file=='na'] <- NA
  cov_file[cov_file=='na'] <- NA
  colnames(pheno_file)[1:3] <- c("FID","IID","count")
  pheno_file[,3] <- as.numeric(pheno_file[,3])
  colnames(cov_file)[1:2] <- c("FID","IID")
  phenodata <- cbind(pheno_file[,1:3],cov_file[,-(1:2)])

  return(list(phenodata=phenodata,genedata=genedata))
}


#' Gene‐based association tests to model zero-inflated count data
#'
#' This function performs gene‐based association tests and omnibus tests
#' between a set of SNPs/genes and zero-inflated count data
#' using ZIP regression or ZINB regression or two-stage SKAT model framework.
#'
#' @param phenodata a data frame containing family and individual IDs for all objects as well as
#' zero-inflated counts as a phenotype and a set of covariates.
#' Each row represents a different individual.
#' The first two columns are Family ID (FID) and Individual ID (IID) respectively.
#' There must be one and only one phenotype in the third column and
#' the phenotype have to be zero-inflated count data which should be non-negative integers, e.g. neuritic plaque counts.
#' Each of the rest of columns represents a different covariate, e.g. age, sex, etc.
#' @param genedata a data frame containing family and individual IDs for all objects as well as numeric genotype data.
#' Each row represents a different individual.
#' The first two columns are Family ID (FID) and Individual ID (IID) respectively.
#' Each of the rest columns represents a seperate gene/SNP marker.
#' The genotype should be coded as 0, 1, 2 and NA for AA, Aa, aa and missing.
#' Both of Family ID (FID) and Individual ID (IID) for each row in the 'genedata'
#' derived from the PLINK formatted files should be in the same order as in the 'phenodata'.
#' The number of rows in 'genedata' should be equal to the number of rows in 'phenodata'.
#' @param name a character string of the name of a gene, e.g. "CETP". The name is case-sensitive.
#' @param weights a character string of pre-specified variant weighting schemes (default="Equal").
#' "Equal" represents no weight,
#' "MadsenBrowning" represents the Madsen and Browning (2009) weight,
#' "Beta" represents the Beta weight.
#' @param missing_cutoff a cutoff of the missing rates of SNPs (default=0.15).
#' Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.
#' @param max_maf a cutoff of the maximum minor allele frequencies (MAF) (default=1, no cutoff).
#' Any SNPs with MAF > cutoff will be excluded from the analysis.
#' @param model character specification of zero-inflated count model family (default="zip").
#' "zip" represents Zero-Inflated Poisson model,
#' "zinb" represents Zero-Inflated Negative Binomial model,
#' "skat" represents the two-stage Sequence Kernel Association Test method.
#' @returns a list of 11 items including the name of gene,
#' the number of rare variants in the gene,
#' the kind of method used for modeling,
#' and individual p-values of gene‐based association tests (burden test and kernel test for both parameters)
#' and omnibus tests using different methods (Cauchy combination test defaulted).
#'
#' \item{GeneName}{the name of gene.}
#' \item{No.Var}{the number of rare variants in the gene.}
#' \item{Method}{the method used to compute the p-values.}
#' \item{p.value_phi_burden / p.value_pi_burden}{single p-value for parameter \eqn{\phi} or \eqn{\pi} using burden test.}
#' \item{p.value_lambda_burden / p.value_mu_burden}{single p-value for parameter \eqn{\lambda} or \eqn{\mu} using burden test.}
#' \item{p.value_phi_kernel / p.value_pi_burden}{single p-value for parameter \eqn{\phi} or \eqn{\pi} using kernel test.}
#' \item{p.value_lambda_kernel / p.value_mu_burden}{single p-value for parameter \eqn{\lambda} or \eqn{\mu} using kernel test.}
#' \item{p.value_burden_combined}{Combined p-value of testing the joint effect of both parameters from burden test using Cauchy combination test (Cauchy-p).}
#' \item{p.value_kernel_combined}{Combined p-value of testing the joint effect of both parameters from kernel test using Cauchy combination test (Cauchy-p).}
#' \item{p.value_overall}{Combined p-value of testing the overall association using Cauchy combination test.}
#'
#' @references Fan, Q., Sun, S., & Li, Y.‐J. (2021). Precisely modeling zero‐inflated count phenotype for rare variants. Genetic Epidemiology, 1–14.
#' @importFrom pscl zeroinfl
#' @importFrom CompQuadForm davies
#' @importFrom stats cov
#' @importFrom stats dbeta
#' @importFrom stats pchisq
#' @importFrom SKAT SKAT_Null_Model
#' @importFrom SKAT SKAT
#' @importFrom RNOmni RankNorm
#' @export

zimfrv <- function(phenodata, genedata, genename = "NA", weights = "Equal", missing_cutoff = 0.15, max_maf = 1, model = "zip", omnibus = "Cauchy"){

  n <- nrow(genedata)
  if(nrow(phenodata) != n){
    stop("phenodata and genedata do not have the same number of rows")
  }
  if((sum(phenodata[,1]==genedata[,1]) != n) | (sum(phenodata[,2]==genedata[,2]) != n)){
    stop("the IDs for each row in phenodata and genedata are not in the same orders")
  }

  phenodata <- phenodata[,-(1:2)] # remove
  genedata <- genedata[,-(1:2)]

  missing_check_ind <- is.na(genedata)
  missing_check <- apply(missing_check_ind, 2, sum)/n
  id_check <- which(missing_check < missing_cutoff)
  genedata_check <- as.matrix(genedata[,id_check])

  sam_maf <- apply(genedata_check, 2, sum, na.rm=T)/(2*n)
  id_rare <- which(sam_maf < max_maf & sam_maf > 0)
  genedata_rare <- as.matrix(genedata_check[,id_rare])
  n_rare <- length(id_rare)
  maf_rare_sam <- sam_maf[id_rare]

  if(weights=="Equal"){
    wt <- rep(1,n_rare)
  }else if(weights=="MadsenBrowning"){
    wt <- 1/sqrt(maf_rare_sam*(1- maf_rare_sam))
  }else if(weights=="Beta"){
    wt <- dbeta(maf_rare_sam,1,25)
  }else{
    stop("invalid weights")
  }
  wt_matrix <- diag(x=wt^2,nrow = n_rare,ncol = n_rare)

  colnames(phenodata)[1] <- "y"

  if(model=="zip"){
    s <- U_fi_lmd(phenodata, genedata_rare)
    s_fi <- s$u_fi
    s_lambda <- s$u_lambda

    s_fi[is.na(s_fi)]=0
    s_lambda[is.na(s_lambda)]=0


    p_fi_burden <- p_burden_single(wt, genedata_rare, s_fi)
    p_lambda_burden <- p_burden_single(wt, genedata_rare, s_lambda)
    p_fi_kernel <- p_kernel_single(wt_matrix, genedata_rare, s_fi)
    p_lambda_kernel <- p_kernel_single(wt_matrix, genedata_rare, s_lambda)

    p_burden_cauchy <- cauchyp(c(p_fi_burden,p_lambda_burden))
    p_kernel_cauchy <- cauchyp(c(p_fi_kernel,p_lambda_kernel))

    p_overall <- cauchyp(c(p_burden_cauchy,p_kernel_cauchy))

    output <- data.frame(GeneName=genename,No.Var=n_rare,Method="ZIP",
                         p_phi_burden=p_fi_burden,p_lambda_burden=p_lambda_burden,
                         p_phi_kernel=p_fi_kernel,p_lambda_kernel=p_lambda_kernel,
                         p_burden_cauchy=p_burden_cauchy,p_kernel_cauchy=p_kernel_cauchy,
                         p_overall=p_overall)


  }else if(model=="zinb"){
    s <- U_phi_mu4zinb(phenodata, genedata_rare)
    s_phi <- s$u_phi
    s_mu <- s$u_mu

    s_phi[is.na(s_phi)]=0
    s_mu[is.na(s_mu)]=0


    p_phi_burden <- p_burden_single(wt, genedata_rare, s_phi)
    p_mu_burden <- p_burden_single(wt, genedata_rare, s_mu)
    p_phi_kernel <- p_kernel_single(wt_matrix, genedata_rare, s_phi)
    p_mu_kernel <- p_kernel_single(wt_matrix, genedata_rare, s_mu)

    p_burden_cauchy <- cauchyp(c(p_phi_burden,p_mu_burden))
    p_kernel_cauchy <- cauchyp(c(p_phi_kernel,p_mu_kernel))

    p_overall <- cauchyp(c(p_burden_cauchy,p_kernel_cauchy))

    output <- data.frame(GeneName=genename,No.Var=n_rare,Method="ZINB",
                         p_phi_burden=p_phi_burden,p_mu_burden=p_mu_burden,
                         p_phi_kernel=p_phi_kernel,p_mu_kernel=p_mu_kernel,
                         p_burden_cauchy=p_burden_cauchy,p_kernel_cauchy=p_kernel_cauchy,
                         p_overall=p_overall)
  }else if(model=="skat"){
    y <- phenodata[,1]
    X <- as.matrix(phenodata[,-1])

    y1 <- ifelse(y>0, 0, 1)
    obj <- SKAT_Null_Model(y1 ~ X, out_type="D")

    id_y_pos <- which(y>0)
    y2 <- y[id_y_pos]
    X2 <- X[id_y_pos,]
    Z <- as.matrix(genedata_rare[id_y_pos,])
    y2t <- RankNorm(y2)
    obj.c <- SKAT_Null_Model(y2t ~ X2, out_type="C")

    p_pi_burden <- SKAT(genedata_rare, obj, weights = wt, r.corr = 1)$p.value
    p_mu_burden <- SKAT(Z, obj.c, weights = wt, r.corr = 1)$p.value
    p_pi_kernel <- SKAT(genedata_rare, obj, weights = wt)$p.value
    p_mu_kernel <- SKAT(Z, obj.c, weights = wt)$p.value

    p_burden_cauchy <- cauchyp(c(p_pi_burden,p_mu_burden))
    p_kernel_cauchy <- cauchyp(c(p_pi_kernel,p_mu_kernel))

    p_overall <- cauchyp(c(p_burden_cauchy,p_kernel_cauchy))

    output <- data.frame(GeneName=genename,No.Var=n_rare,Method="SKAT",
                         p_pi_burden=p_pi_burden,p_mu_burden=p_mu_burden,
                         p_pi_kernel=p_pi_kernel,p_mu_kernel=p_mu_kernel,
                         p_burden_cauchy=p_burden_cauchy,p_kernel_cauchy=p_kernel_cauchy,
                         p_overall=p_overall)

  }else{
    stop("invalid model family")
  }

  return(output)

}
