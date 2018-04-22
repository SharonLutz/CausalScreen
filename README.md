## CausalScreen
The CausalScreen R package implements a screening step to increase power when testing for direct genetic effects of multiple SNPs in family based association studies using causal inference methodology. 

#### Installation
```
install.packages("devtools") #The devtools package must be installed first

devtools::install_github("SharonLutz/CausalScreen")
```
#### Example
For the sample dataset dataS, adjusting for age and gender and given the parental genotypes (i.e. P1 and P2), one can test if each  of the 100 SNPs (i.e. X) is directly associated with FEV given that the SNP is associated with FEV through the intermediate variable height and weight. The code to load the sample dataset and run CausalScreen are given below.
```
library(CausalScreen)
?CausalScreen # For details on this function and how to choose input variables

data("dataS")
X<-dataS[,1:100]
P1<-dataS[,101:200]
P2<-dataS[,201:300]
height<-dataS[,"height"]
weight<-dataS[,"weight"]
fev<-dataS[,"fev"]
age<-dataS[,"age"]
gender<-dataS[,"gender"]
CausalScreen(X,P1,P2,height,weight,fev,cbind(age,gender))
```

#### Output
For this example, a subset of the output is given below. The first 2 SNPs directly act on FEV given the alternative pathways through height and weight after adjusting for age and gender (p-value=2.7E-8 and 2.7E-6). For more details on this method, please see the references listed below.

```
     snpOrder ScreeningStepOrder   ScreenStat    Statistic       pvalue
x1          1                  1 2.081983e+00 3.093912e+01 2.662510e-08
x2          2                  2 5.603640e-01 2.204155e+01 2.668119e-06
x78        78                  3 1.137398e-01 1.543370e-01 6.944247e-01
x67        67                  4 1.101877e-01 1.518423e+00 2.178583e-01
x40        40                  5 8.921338e-02 5.541485e-01 4.566276e-01
.
.
.
```
#### References
**Lutz SM**, Vansteelandt S, Lange C. (2013) Testing for Direct Genetic Effects Using a Screening Step in Family-Based Association Studies. *Frontiers in Genetics*. 4 (243).

Vansteelandt S, Goetgeluk S, **Lutz S**, Waldamn I, Lyon H, Schadt EE, Weiss ST, Lange C. (2009) On the Adjustment for Covariates in Genetic Association Analysis: A Novel, Simple Principle to Infer Direct Effects. *Genetic Epidemiology*. 33(5): 394-405.
