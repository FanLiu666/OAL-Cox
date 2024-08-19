## Writen for R version 4.3.3
library(MASS) # version 3.3.1
library(SurvMetrics)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(openxlsx)
library(OALCox)

source('AKM_rmst.R')
source('sim_data.R')
source('AL-Cox.R')

############################ Compare RASCE estimated by differet methods ############################
#' RASCE_comparison function
#'
#' @param rho: correlation of covariates
#' @param n: sample size
#' @param p: number of covariates
#' @param theta: Censored rate
#' @param tau: time point of interest
#' @param rep_num: number of repetitions
#' @param lambda_vec: regularization parameter for OAL-Cox model
#' @return RASCE_comparison plot and the model evaluation index

RASCE_comparison <- function(rho, n, p, theta_tau_pairs, rep_num, lambda_vec, compare_table) {
  for (pair in theta_tau_pairs) {
    theta <- pair[[1]]
    tau <- pair[[2]]
    for (r in rho) {
      for(t in p) {
        ### calculate the true value of RASCE
        RASCE_original <- c()
        for(i in 1:500){
          a <- simulate_data(mean_x = 0,
                             sig_x = 1,
                             rho = r,
                             n = 10000,
                             p = t,
                             lambda = 0.0001,
                             alpha = 3,
                             theta = theta)
          Data <- a[[1]]
          var.list <- a[[2]]

          # calculate
          RASCE_True <- calc_ipw_rmst(data = Data,
                                      A = Data$A,
                                      time = Data$time,
                                      status = Data$status,
                                      covariates = var.list,
                                      tau = tau)
          RASCE_original[i] <- RASCE_True
        }

        True_RASCE <- round(mean(RASCE_original), 3)

        for (s in n) {
          res_df <- data.frame(RASCE = numeric(), method = character(), circle = numeric())

          for (i in 1:rep_num) {
            ### simulate data
            a <- simulate_data(mean_x = 0,
                               sig_x = 1,
                               rho = r,
                               n = s,
                               p = t,
                               lambda = 0.0001,
                               alpha = 3,
                               theta = theta)
            Data <- a[[1]]
            var.list <- a[[2]]

            ### run OAL-Cox model
            b <- OAL_Cox(data = Data,
                         A = Data$A,
                         time = Data$time,
                         status = Data$status,
                         covariates = var.list,
                         lambda_vec = lambda_vec,
                         gamma_convergence_factor = 5,
                         tau = tau)
            # RASCE for OAL-Cox model
            RASCE_OAL_Cox <- b[[1]]

            ### run AL-Cox model
            d <- AL_Cox(data = Data,
                        A = Data$A,
                        time = Data$time,
                        status = Data$status,
                        covariates = var.list,
                        lambda_vec = lambda_vec,
                        tau = tau)
            RASCE_AL_Cox <- d[[1]]

            ### run conf weighting model
            var_Conf <- grepl("Xc", names(Data))
            var_Conf <- names(Data)[var_Conf]
            RASCE_Conf <- calc_ipw_rmst(data = Data,
                                        A = Data$A,
                                        time = Data$time,
                                        status = Data$status,
                                        covariates = var_Conf,
                                        tau = tau)

            ### run Targ weighting model
            var_Targ <- grepl("Xc", names(Data)) | grepl("Xp", names(Data))
            var_Targ <- names(Data)[var_Targ]
            RASCE_Targ <- calc_ipw_rmst(data = Data,
                                        A = Data$A,
                                        time = Data$time,
                                        status = Data$status,
                                        covariates = var_Targ,
                                        tau = tau)

            ### run PreT weighting model
            var_PreT <- grepl("Xc", names(Data)) | grepl("Xi", names(Data))
            var_PreT <- names(Data)[var_PreT]
            RASCE_PreT <- calc_ipw_rmst(data = Data,
                                        A = Data$A,
                                        time = Data$time,
                                        status = Data$status,
                                        covariates = var_PreT,
                                        tau = tau)

            ### run PotConf weighting model
            var_PotConf <- grepl("Xc", names(Data)) | grepl("Xp", names(Data)) | grepl("Xi", names(Data))
            var_PotConf <- names(Data)[var_PotConf]
            RASCE_PotConf <-calc_ipw_rmst(data = Data,
                                          A = Data$A,
                                          time = Data$time,
                                          status = Data$status,
                                          covariates = var_PotConf,
                                          tau = tau)

            # combind RASCE for different method
            res_df <- rbind(res_df, data.frame(RASCE = c(RASCE_OAL_Cox,
                                                         RASCE_AL_Cox,
                                                         RASCE_Conf,
                                                         RASCE_Targ,
                                                         RASCE_PreT,
                                                         RASCE_PotConf),
                                               method = c("OAL-Cox", "AL-Cox", "Conf", "Targ", "PreT", "PotConf"),
                                               circle = rep(i, 6)))
          }

          # specifies the abscissa order
          model <- c("OAL-Cox", "AL-Cox", "Targ", "Conf", "PreT", "PotConf")
          res_df$method <- factor(res_df$method, levels = model)

          ### plot the RASCE_compare figure
          plot <- ggplot(res_df, aes(x = method, y = RASCE))+
            stat_boxplot(geom = "errorbar", width = 0.1, size = 0.8)+
            geom_boxplot()+
            xlab("")+
            ylab("Estimate")+
            geom_hline(yintercept = True_RASCE, linetype = 'dotted', col = 'red')+
            theme_bw()+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

          file_name <- paste("~/RASCE_method1/", "1_",
                             s, "_",
                             t, "_",
                             r, "_",
                             theta, ".png", sep = "")
          ggsave(file_name, plot)

          ### save the model evaluation index
          i <- 1
          RB <- c()
          MSE <- c()
          for (meth in model) {
            # calculate relative risk
            RASCE <- res_df[res_df$method == meth, 1]
            estimate <- mean(res_df[res_df$method == meth, 1])
            RB_i <- (estimate-True_RASCE)/True_RASCE
            RB[i] <- RB_i

            # calculate mean squared error
            MSE_i <- sum((RASCE - True_RASCE)^2)/rep_num
            MSE[i] <- MSE_i

            i <- i+1
          }

          name <- paste("_n_", s,
                        "_p_", t,
                        "_rho_", r,
                        sep = "")
          compare_table[compare_table$cencored_rate == theta, name] <- c(RB, MSE)

        }
      }
    }
  }
  write.xlsx(compare_table, paste("~/RASCE_method1/", "n_", n, ".xlsx", sep = ""))
}

# simulate
model <- c("OAL-Cox", "AL-Cox", "Targ", "Conf", "PreT", "PotConf")
lambda_vec <- c(-10,-5,-2,-1,-0.75,-0.5,-0.25, 0.25, 0.49)
compare_index <- data.frame(cencored_rate = c(rep(50,12), rep(100,12)),
                            evaluation_type = c(rep('RB', 6), rep('MSE', 6), rep('RB', 6), rep('MSE', 6)),
                            method = c(rep(model, 4)))

rho_values <- c(0.2, 0.5)
theta_tau_pairs <- list(c(50, 30), c(100, 60))
p <- c(100, 200)
n <- c(500, 1000)
rep <- 500
RASCE_comparison(rho = rho_values,
                 n = n,
                 p = p,
                 theta_tau_pairs =theta_tau_pairs,
                 rep_num = rep,
                 lambda_vec = lambda_vec,
                 compare_table = compare_index)
