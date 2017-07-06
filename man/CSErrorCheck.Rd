\name{CSErrorCheck}
\alias{CSErrorCheck}

\title{CSErrorCheck}
\description{
Performs error checks on the input for the CausalScreen function.
}

\usage{
CSErrorCheck(X, P1, P2, l, k, y, z = NULL)
}
\arguments{
  \item{X}{
  is a matrix of the child's genotype.
}
  \item{P1}{
  is a matrix of parent 1's genotype (ie; the father).
}
  \item{P2}{
  is a matrix of parent 2's genotype (ie; the mother).
}
  \item{l}{
  is an intermediate phenotype such as height.
}
  \item{k}{
  is the intermediate phenotype of interest such as weight.
}
  \item{y}{
  is the normally distributed outcome of interest such as FEV.
}
  \item{z}{
  is a matrix of covariates to adjust for, such as age and gender.
}
}

\value{
This function exits with an error message if the error checks are not met.
}

\author{
Sharon Lutz, Annie Thwing
}
