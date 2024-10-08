
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ZIM4rv.R
\name{zimfrv}
\alias{zimfrv}
\title{Gene‐based association tests to model zero-inflated count data}
\usage{
zimfrv(
  phenodata,
  genedata,
  name = NULL,
  weights = "Equal",
  missing_cutoff = 0.15,
  max_maf = 1,
  model = "zip"
)
}
\arguments{
\item{phenodata}{a data frame containing family and individual IDs for all objects as well as
zero-inflated counts as a phenotype and a set of covariates.
Each row represents a different individual.
The first two columns are Family ID (FID) and Individual ID (IID) respectively.
There must be one and only one phenotype in the third column and
the phenotype have to be zero-inflated count data which should be non-negative integers, e.g. neuritic plaque counts.
Each of the rest of columns represents a different covariate, e.g. age, sex, etc.}

\item{genedata}{a data frame containing family and individual IDs for all objects as well as numeric genotype data.
Each row represents a different individual.
The first two columns are Family ID (FID) and Individual ID (IID) respectively.
Each of the rest columns represents a seperate gene/SNP marker.
The genotype should be coded as 0, 1, 2 and NA for AA, Aa, aa and missing.
Both of Family ID (FID) and Individual ID (IID) for each row in the 'genedata'
derived from the PLINK formatted files should be in the same order as in the 'phenodata'.
The number of rows in 'genedata' should be equal to the number of rows in 'phenodata'.}

\item{name}{a character string of the name of a gene, e.g. "CETP". The name is case-sensitive.}

\item{weights}{a character string of pre-specified variant weighting schemes (default="Equal").
"Equal" represents no weight,
"MadsenBrowning" represents the Madsen and Browning (2009) weight,
"Beta" represents the Beta weight.}

\item{missing_cutoff}{a cutoff of the missing rates of SNPs (default=0.15).
Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.}

\item{max_maf}{a cutoff of the maximum minor allele frequencies (MAF) (default=1, no cutoff).
Any SNPs with MAF > cutoff will be excluded from the analysis.}

\item{model}{character specification of zero-inflated count model family (default="zip").
"zip" represents Zero-Inflated Poisson model,
"zinb" represents Zero-Inflated Negative Binomial model,
"skat" represents the two-stage Sequence Kernel Association Test method.}
}
\value{
a list of 10 items including the name of gene, the number of rare variants in the gene
and p-values of gene‐based association tests (burden test and kernel test for both parameters)
and omnibus tests using different methods (Cauchy’s and Fisher’s respectively).

\item{GeneName}{the name of gene.}
\item{No.Var}{the number of rare variants in the gene.}
\item{p.value_phi_burden / p.value_pi_burden}{single p-value for parameter \eqn{\phi} or \eqn{\pi} using burden test.}
\item{p.value_lambda_burden / p.value_mu_burden}{single p-value for parameter \eqn{\lambda} or \eqn{\mu} using burden test.}
\item{p.value_phi_kernel / p.value_pi_burden}{single p-value for parameter \eqn{\phi} or \eqn{\pi} using kernel test.}
\item{p.value_lambda_kernel / p.value_mu_burden}{single p-value for parameter \eqn{\lambda} or \eqn{\mu} using kernel test.}
\item{p.value_burden_cauchy}{Combined p-value of testing the joint effect of both parameters from burden test using Cauchy combination test (Cauchy-p).}
\item{p.value_burden_fisher}{Combined p-value of testing the joint effect of both parameters from burden test using Fisher's method (Fisher-p).}
\item{p.value_kernel_cauchy}{Combined p-value of testing the joint effect of both parameters from kernel test using Cauchy combination test (Cauchy-p).}
\item{p.value_kernel_fisher}{Combined p-value of testing the joint effect of both parameters from kernel test using Fisher's method (Fisher-p).}
}
\description{
This function performs gene‐based association tests and omnibus tests
between a set of SNPs/genes and zero-inflated count data
using ZIP regression or ZINB regression or two-stage SKAT model framework.
}
\references{
Fan, Q., Sun, S., & Li, Y.‐J. (2021). Precisely modeling zero‐inflated count phenotype for rare variants. Genetic Epidemiology, 1–14.
}
