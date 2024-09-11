## Writen for R version 4.3.3
library(MASS) # version 3.3.1
library(SurvMetrics)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(openxlsx)
library(OALCox)

source('sim_data.R')
source('AL-Cox.R')

############################ percentage of covariance selected ############################
#' variable_select function
#'
#' @param rho: correlation of covariates
#' @param n: sample size
#' @param p: number of covariates
#' @param theta: Censored rate
#' @param tau: time point of interest
#' @param lambda_vec: regularization parameter for OAL-Cox model
#' @param rep_num: number of repetitions
#' @return variable selected proportion plot

variable_select <- function(rho, n_p_pairs, theta_tau_pairs, lambda_vec, rep_num) {
  for (pair in theta_tau_pairs) {
    theta <- pair[[1]]
    tau <- pair[[2]]

    for(pair2 in n_p_pairs) {
      n <- pair2[[1]]
      p <- pair2[[2]]

      for (r in rho) {
        var_select <- data.frame(Proportion = numeric(), index = numeric(), Group = character())

        covar_select_OAL_Cox <- data.frame(matrix(nrow = rep_num, ncol = p))
        covar_select_AL_Cox <- data.frame(matrix(nrow = rep_num, ncol = p))

        for (i in 1:rep_num) {
          a <- simulate_data(mean_x = 0,
                             sig_x = 1,
                             rho = r,
                             n = n,
                             p = p,
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
                       gamma_convergence_factor = 10,
                       tau = tau)
          coeff_OAL_Cox <- b[[3]]
          covar_select_OAL_Cox[i,] <- t(coeff_OAL_Cox)

          ### run AL-Cox model
          d <- AL_Cox(data = Data,
                      A = Data$A,
                      time = Data$time,
                      status = Data$status,
                      covariates = var.list,
                      lambda_vec = lambda_vec,
                      tau = tau)
          coeff_AL_Cox <- d[[2]]
          covar_select_AL_Cox[i,] <- t(coeff_AL_Cox)
        }

        # calculate the variable-seletion proportion of OAL-Cox
        c1 <- apply(covar_select_OAL_Cox, 2, sum) / rep_num
        var_select <- rbind(var_select, data.frame(Proportion = c1[1:20],
                                                   index = c(1:20),
                                                   Group = rep("OAL-Cox", 20)))

        # calculate the variable-seletion proportion of AL-Cox
        c2 <- apply(covar_select_AL_Cox, 2, sum) / rep_num
        var_select <- rbind(var_select, data.frame(Proportion = c2[1:20],
                                                   index = c(1:20),
                                                   Group = rep("AL-Cox", 20)))

        # calculate the variable-seletion proportion of reference method
        c3 <- c(rep(1,4), rep(0, p-4))
        var_select <- rbind(var_select, data.frame(Proportion = c3[1:20],
                                                   index = c(1:20),
                                                   Group = rep("Reference", 20)))

        # plot the variable-seletion proportion of different method
        var_select$Group <- factor(var_select$Group, levels = c("OAL-Cox", "AL-Cox", "Reference"))

        plot <- ggplot(data = var_select, aes(x = index, y = Proportion, linetype = Group, color = Group))+
          geom_line(size = 0.8)+
          xlab("Covariate index")+
          ylab("Proportion of times covariate selected")+
          geom_hline(yintercept = 0.9, linetype = 'dotted', col = 'black')+
          theme_bw()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = c(1, 1),
                legend.justification = c(1, 1),
                legend.box.background = element_rect(color="black"))+
          scale_x_continuous(limits = c(0,20), breaks = seq(0, 20, 5))

        file_name <- paste("~/variable_select_1/", "1_",
                           n, "_",
                           r, "_",
                           theta, ".png", sep = "")
        ggsave(file_name, plot)

      }
    }
  }
}

# simulate
lambda_vec <- c(-10,-5,-2,-1,-0.75,-0.5,-0.25, 0.25, 0.49)
rho_values <- c(0, 0.2)
n_p_pairs <- list(c(500, 50), c(1000, 50))
theta_tau_pairs <- list(c(50, 30), c(100, 60))
rep_num <- 500
variable_select(rho = rho_values,
                n_p_pairs = n_p_pairs,
                theta_tau_pairs =theta_tau_pairs,
                lambda_vec = lambda_vec,
                rep_num = rep_num)
