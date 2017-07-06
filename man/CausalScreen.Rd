\name{CausalScreen}
\alias{CausalScreen}
\title{CausalScreen}
\description{
CausalScreen implements a screening step to increase power when testing for direct genetic effects of multiple SNPs in family based association studies using causal inference methodology.
}
\usage{
CausalScreen(X, P1, P2, l, k, y, z = NULL)
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
\item{snpOrder}{The order the SNP appears in the matrix}
\item{ScreeningStepOrder}{The order of the SNP for the screening step.}
\item{ScreenStat}{The statistic from the screeningstat that is not based on the child's genotype.}
\item{Statistic and p-value}{Tests the direct genetic effect of the genotype of the SNP x on the outcome of interest y, given the intermediate phenotypes l and k.}
}

\references{
Lutz SM, Vansteelandt S, Lange C. (2013) Testing for Direct Genetic Effects Using a Screening Step in Family-Based Association Studies. Frontiers in Genetics. 4 (243).

Vansteelandt S, Goetgeluk S, Lutz S, Waldamn I, Lyon H, Schadt EE, Weiss ST, Lange C. (2009) On the Adjustment for Covariates in Genetic Association Analysis: A Novel, Simple Principle to Infer Direct Effects. Genetic Epidemiology. 33(5): 394-405.
}

\author{
Sharon Lutz, Annie Thwing
}

\keyword{ causal }
\keyword{ screen }
\keyword{ SNP }