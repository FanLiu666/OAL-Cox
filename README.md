# OAL-Cox: Outcome-adaptive Lasso for Cox Proportional Hazard Model
The OAL-Cox model was designed for accurate variable selection and casual estimation in survival data. Confounders and outcome predictors are selected by minimizing weighted absolute mean difference (wAMD) considereing covariates balance. Simultaneously, we utilize coefficient estimates from the
minimizer of the negative log Cox partial likelihood function to construct penalty weights, and employ the restricted average survival time causal effect (RASCE) for causal effect estimation in survival data.
## Installation
You can install the development version of OAL-Cox from
[GitHub](https://github.com/) with:

``` r
library(devtools)
install_local(“/path/to/lqa_1.0-3.tar.gz”)
install_github("FanLiu666/OAL-Cox")
```
## Example   

``` r
library(OALCox)
library(joint.Cox)
#> Loading required package: lqa

# import data
data("dataOvarian1")
data <- dataOvarian1
lambda_vec <- c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
gamma <- 5

# data processing
Data_use <- data[, -3]
colnames(Data_use)[1:3] <- c("time", "status", "A")
Data <- Data_use[Data_use$time != 0, ]
var.list <- colnames(Data_use[, -c(1:3)])

OAL_Cox(data = Data,
        A = Data$A, 
        time = Data$time, 
        status = Data$status, 
        covariates = var.list, 
        lambda_vec = lambda_vec, 
        gamma_convergence_factor = gamma,
        tau = 1095)
```
