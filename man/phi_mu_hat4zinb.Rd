
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ZIM4rv.R
\name{phi_mu_hat4zinb}
\alias{phi_mu_hat4zinb}
\title{Estimation of phi_hat, mu_hat and alpha_hat for ZINB model}
\usage{
phi_mu_hat4zinb(simud)
}
\arguments{
\item{simud}{a data frame containing a phenotype named y and covariates}
}
\value{
a list of 3 estimations of parameters for each subject
}
\description{
This function gives the estimation of 3 parameters phi, mu and alpha in ZINB model
for each subject under the null hypothesis.
}
\details{
This function first fits zero‐inflated negative binomial regression of phenotype y
on the covariates only to obtain the estimates of regression coefficients
and inverse dispersion
and then compute the estimations of phi, mu and alpha.
}
\seealso{
\code{zeroinfl}
}
