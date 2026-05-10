# Bayesian Hierarchical Analysis: County-Level COVID Vaccination Rates

Fit a nested Bayesian Beta-Binomial hierarchical model (counties within states) 
to a cross-sectional snapshot of COVID-19 vaccination completion counts across 
2,921 U.S. counties — borrowing strength across counties while preserving 
local heterogeneity.

**Why Bayesian:** County-level vaccination data has highly variable denominators 
(small counties produce noisy proportions). Partial pooling stabilizes estimates 
for small counties while letting large counties be driven by local data. 
OLS on raw proportions would produce unreliable estimates for low-population counties.

**Key finding:** Hierarchical model outperformed independent Beta(1,1) baseline 
by 3,292 WAIC points. Shrinkage plots confirm expected partial pooling behavior.

**Tools:** R (tidyverse, ggplot2, matrixStats)  
**Methods:** Bayesian hierarchical modeling, Beta-Binomial, grid approximation, 
WAIC model comparison, posterior predictive checking, partial pooling and when to use it  
**Course:** STP 505 — Bayesian Statistics, ASU (Spring 2026)
