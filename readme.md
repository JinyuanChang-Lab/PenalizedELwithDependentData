## Official implementations of `Empirical likelihood approach for high-dimensional
moment restrictions with dependent data`

## Introduction

We introduce an empirical likelihood approach for high-dimensional moment restrictions with dependent data. Since the numbers of parameters and moment conditions can potentially exceed the sample size, we add a double penalty to the empirical likelihood criterion function to introduce sparsity for dimensional reduction.  Under certain regularity conditions, we establish asymptotic guarantees of our proposed method, thereby making it an attractive option for estimating and inferring large multivariate time series models, including VAR, LP, and MARCH/MGARCH. The versatility of our method is illustrated through extensive Monte Carlo simulations and three empirical applications, which include analyzing US sectoral inflation rates, fiscal policy multipliers, and the volatility correlation within China’s banking sector.

## Folders

* simulation: includes codes and data for Section 5 （Simulations）in the paper.

* realdata: includes codes and data for Section 6（Empirical applications）in the paper.

## Codes

The codes in the folders include the implementation of the proposed estimation and inference methods, called PEL and PPEL in the paper, respectively, and their applications to VAR, LP, and MGARCH models.