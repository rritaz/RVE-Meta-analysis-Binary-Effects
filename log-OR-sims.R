# Load libraries
library(tidyverse)
library(readxl)

source("df_HC2.R")
source("df_HC3.R")

# Read in excel file containing n values
n_values <- read_excel("n_values.xlsx", col_names = FALSE)

pi_c <- c(0.06, 0.1, 0.5, 0.3) # True rate in control group
pi_t <- c(0.08, 0.1, 0.5, 0.4) # True rate in treatment group

# Values of constants A and B in Hartung RE procedure
A <- 0.95
B <- 1.05

#----------------------------------------------------------------------------#
# Data generation process: 2 steps, study specific pi's and observed p's

## Step 1
get_pi_tk <- function(alpha_t, beta_t){
  return(rbeta(1, alpha_t, beta_t))
}

get_pi_ck <- function(alpha_c, beta_c){
  return(rbeta(1, alpha_c, beta_c))
}

## step 2
get_m_tk <- function(n, pi_tks){
  return(rbinom(1, as.numeric(n), as.numeric(pi_tks)))
}

get_m_ck <- function(n, pi_cks){
  return(rbinom(1, as.numeric(n), as.numeric(pi_cks)))
}

# Generating data
generate_pi_tk_forCol <- function(kVal, alpha_t, beta_t){
  pi_tks <- vector()
  for(i in 1:kVal) {
    pi_tks[i] <- get_pi_tk(alpha_t, beta_t)
  }
  return(pi_tks)
}

generate_pi_ck_forCol <- function(kVal, alpha_c, beta_c){
  pi_cks <- vector()
  for(i in 1:kVal) {
    pi_cks[i] <- get_pi_ck(alpha_c, beta_c)
  }
  return(pi_cks)
}

generate_m_tk_forCol <- function(colNum, kVal, pi_tks){
  m_tks <- vector()
  for(i in 1:kVal) {
    m_tks[i] <- get_m_tk(n_values[i, colNum], pi_tks[i])
  }
  return(m_tks)
}

generate_m_ck_forCol <- function(colNum, kVal, pi_cks){
  m_cks <- vector()
  for(i in 1:kVal) {
    m_cks[i] <- get_m_ck(n_values[i, colNum], pi_cks[i])
  }
  return(m_cks)
}

log_odds_ratioforCol <- function(kVal, colNum, m_tks, m_cks){
  log_odds_ratios <- vector()
  for(i in 1:kVal){
    log_odds_ratios[i] <- log(((m_tks[i]+0.5)/((as.numeric(n_values[i, colNum])-m_tks[i])+0.5))*
                                (((as.numeric(n_values[i, colNum])-m_cks[i])+0.5)/(m_cks[i]+0.5)))
  }
  return(log_odds_ratios)
} # Log Odds Ratio for k studies

generate_sigma_forCol <- function(kVal, colNum, m_tks, m_cks){
  sigmas <- vector()
  for(i in 1:kVal) {
    sigmas[i] <- (1/(m_tks[i]+0.5)) + (1/((as.numeric(n_values[i, colNum])-m_tks[i])+0.5)) +
      (1/((as.numeric(n_values[i, colNum])-m_cks[i])+0.5)) + (1/(m_cks[i]+0.5))
  }
  return(sigmas) # Variance of log OR assuming equal treatment and control group obs
}
#----------------------------------------------------------------------------#

# Meta-analysis 

weight <- function(sigma) {
  output <- 1 / sigma
  return(output)
} # Fixed effect weights

weightTrue <- function(tau_sqVal, sigma) {
  output <- 1 / (sigma + tau_sqVal)
  return(output)
} # Random effect weights

sumW <- function(kVal, powToRaise, sigmas) {
  sum_w <- 0
  for (i in 1:kVal) {
    sum_w <- sum_w + weight(sigmas[i])^powToRaise
  }
  return(as.numeric(sum_w))
} # Sum of fixed effect weights

sumTrueW <- function(kVal, tauHat_value, sigmas) {
  sum_w <- 0
  for (i in 1:kVal) {
    sum_w <- sum_w + weightTrue(tauHat_value, sigmas[i])
  }
  return(as.numeric(sum_w))
} # Sum of random effect weights

tau_sq_fun <- function(pi_t, pi_c) {
  total <- (pi_t*(1-pi_t)+pi_c*(1-pi_c))
  return(as.numeric(total)/2)
} # True heterogeneity between studies

average_nValues <- function(colNum, kVal) {
  total <- 0
  for(i in 1:kVal) {
    total <- total + (n_values[i, colNum])
  }
  average <- total/kVal
  return(as.numeric(average))
} # Average sample size of k studies

theta_hat <- function(kVal, tauHat_value, log_odds_ratios, sigmas) {
  sum_Numerator <- 0
  for (i in 1:kVal) {
    sum_Numerator <- sum_Numerator + weightTrue(tauHat_value, sigmas[i]) * log_odds_ratios[i]
  }
  return (as.numeric(sum_Numerator/sumTrueW(kVal, tauHat_value, sigmas)))
} # Summary estimate in RE procedure

#----------------------------------------------------------------------------#

Q <- function(log_odds_ratios, thetaHat_FE_val, kVal, sigmas){
  numerator <- 0
  for (i in 1:length(log_odds_ratios)){
    numerator <- numerator + weight(sigmas[i])*(log_odds_ratios[i]-thetaHat_FE_val)^2
  }
  return(numerator)
} # Weighted sum of squares Q using FE weights

theta_hatFE <- function(kVal, log_odds_ratios, sigmas) {
  sum_Numerator <- 0
  for (i in 1:kVal) {
    sum_Numerator <- sum_Numerator + weight(sigmas[i]) * log_odds_ratios[i]
  }
  return (as.numeric(sum_Numerator/sumW(kVal, 1, sigmas)))
} # theta hat with fixed effect weights


denomConst <- function(kVal, sigmas) {
  toReturn <- (sumW(kVal, 1, sigmas) - sumW(kVal, 2, sigmas)/sumW(kVal, 1, sigmas))
  return(as.numeric(toReturn))
} # Constant function of weights in DL denominator of tau hat squared


tauHat_DL <- function(log_odds_ratios, thetaHat_FE_val, kVal, sigmas){
  toReturn <- ((Q(log_odds_ratios, thetaHat_FE_val, kVal, sigmas)-(kVal-1))/denomConst(kVal, sigmas))
  if (toReturn<0){
    return(0)
  }else{
    return(as.numeric(toReturn))
  }
} # DerSimonian and Laird tau hat

betas <- function(kVal, tauHat_value, sigmas){
  beta <- vector()
  for(i in 1:kVal) {
    beta[i] <- weightTrue(tauHat_value, sigmas[i])/sumTrueW(kVal, tauHat_value, sigmas)
  }
  return(as.numeric(beta))
} # Normalized weights

sumOfBetasSQ <- function(betas, kVal){
  sum_betas <- 0
  for(i in 1:kVal){
    sum_betas <- sum_betas + betas[i]^2
  }
  return(as.numeric(sum_betas))
} # Sum of normalized weights

ksi_i_Beta <- function(betas, kVal){
  ksi_i_Beta <- vector()
  for(i in 1:kVal){
    ksi_i_Beta[i] <- betas[i]+((-betas[i]+(betas[i])^2)/(1-sumOfBetasSQ(betas, kVal)))
  }
  return(as.numeric(ksi_i_Beta))
}

QofBeta_term2 <- function(kVal, sigmas, ksi_i_Beta){
  term2 <- 0
  for(i in 1:kVal){
    term2 <- term2 + ksi_i_Beta[i]*sigmas[i]
  }
  return(as.numeric(term2))
}

Q_RE <- function(log_odds_ratios, thetaHat_val, kVal, betas){
  numerator <- 0
  for (i in 1:length(log_odds_ratios)){
    numerator <- numerator + betas[i]*(log_odds_ratios[i]-thetaHat_val)^2
  }
  return(numerator)
} # Weighted sum of squares - using normalized weights

RofBeta <- function(betas, sigmas){
  numerator <- 0
  for (i in 1:length(sigmas)){
    numerator <- numerator + betas[i]^2*sigmas[i]
  }
  return(as.numeric(numerator))
} # Lower estimator

hc1 <- function(log_odds_ratios, thetaHat_val, kVal, sigmas, tauHat_value){
  numerator <- 0
  for (i in 1:length(log_odds_ratios)){
    numerator <- numerator + weightTrue(tauHat_value, sigmas[i])^2*(kVal/(kVal-1))*
      (log_odds_ratios[i]-thetaHat_val)^2
  }
  return(numerator/sumTrueW(kVal, tauHat_value, sigmas)^2)
} 

hc2 <- function(log_odds_ratios, thetaHat_val, kVal, sigmas, tauHat_value){
  numerator <- 0
  for (i in 1:length(log_odds_ratios)){
    numerator <- numerator + weightTrue(tauHat_value, sigmas[i])^2*(1/(1-(weightTrue(tauHat_value, sigmas[i])/
                                                                            sumTrueW(kVal, tauHat_value, sigmas))))*
      (log_odds_ratios[i]-thetaHat_val)^2
  }
  return(numerator/sumTrueW(kVal, tauHat_value, sigmas)^2)
} 

hc3 <- function(log_odds_ratios, thetaHat_val, kVal, sigmas, tauHat_value){
  numerator <- 0
  for (i in 1:length(log_odds_ratios)){
    numerator <- numerator + weightTrue(tauHat_value, sigmas[i])^2*(1/(1-(weightTrue(tauHat_value, sigmas[i])/
                                                                            sumTrueW(kVal, tauHat_value, sigmas))))^2*
      (log_odds_ratios[i]-thetaHat_val)^2
  }
  return(numerator/sumTrueW(kVal, tauHat_value, sigmas)^2)
} 

c_hc2 <- function(kVal, betas){
  c <- vector()
  for (i in 1:kVal){
    c[i] <- (betas[i]^2)/(1-betas[i])
  }
  return(c)
} # c for HC2

c_hc3 <- function(kVal, betas){
  c <- vector()
  for (i in 1:kVal){
    c[i] <- (betas[i]^2)/((1-betas[i])^2)
  }
  return(c)
} # c for HC3

main <- function(){
  kVal <- 10 # Number of studies
  totalNumCol <- length(n_values[1,])

  mainResults <- matrix(NA, nrow=3*totalNumCol*4, ncol=32)
  currentRow <- 1
  for(i in 1:4){ # i is the index for true rate vector (i.e., 4)
    for(j in 1:totalNumCol) {
      phi <- c(0.000001, 1/(1+2*average_nValues(j, kVal)), 1/(1+average_nValues(j, kVal)))
      for(k in 1:length(phi)){
        alpha_t <- pi_t[i]*((1/phi[k])-1) 
        beta_t <- (1-pi_t[i])*((1/phi[k])-1)
        
        alpha_c <-  pi_c[i]*((1/phi[k])-1)
        beta_c <- (1-pi_c[i])*((1/phi[k])-1)
        
        theta <- log((pi_t[i]/(1-pi_t[i]))/(pi_c[i]/(1-pi_c[i]))) # True log OR
        currentTau <- phi[k]*tau_sq_fun(pi_t[i], pi_c[i]) 
        
        pi_tks <- generate_pi_tk_forCol(kVal, alpha_t, beta_t)
        pi_cks <- generate_pi_ck_forCol(kVal, alpha_c, beta_c)
        
        m_tks <- generate_m_tk_forCol(j, kVal, pi_tks)
        m_cks <- generate_m_ck_forCol(j, kVal, pi_cks)
  
        log_odds_ratios <- log_odds_ratioforCol(kVal, j, m_tks, m_cks)
        sigmas <- generate_sigma_forCol(kVal, j, m_tks, m_cks)
  
        thetaHat_FE_val <- theta_hatFE(kVal, log_odds_ratios, sigmas) # theta hat value
        tauHat_value <- tauHat_DL(log_odds_ratios, thetaHat_FE_val, kVal, sigmas)
  
        thetaHat_val <- theta_hat(kVal, tauHat_value, log_odds_ratios, sigmas) # Random effect meta analytic estimate
        
        betas <- betas(kVal, tauHat_value, sigmas)
        c_hc2 <- c_hc2(kVal, betas)
        c_hc3 <- c_hc3(kVal, betas)
        sumOfBetasSQ <- sumOfBetasSQ(betas, kVal)
        lambdaOfBeta <- sumOfBetasSQ(betas, kVal)/(1-sumOfBetasSQ(betas, kVal))
        ksi_i_Beta <- ksi_i_Beta(betas, kVal)
  
        QofBeta <- lambdaOfBeta*Q_RE(log_odds_ratios, thetaHat_val, kVal, betas) + QofBeta_term2(kVal, sigmas, ksi_i_Beta)
        RofBeta <- RofBeta(betas, sigmas)
        LofBeta <- min(1, max(0, ((QofBeta/RofBeta) - A)/(B-A)))
        qofBeta <- LofBeta*QofBeta + (1-LofBeta)*RofBeta # Truncated Hartung variance estimator
        
        Var_QofBeta <- (lambdaOfBeta^2)*(1/sumTrueW(kVal, tauHat_value, sigmas))^2*(2*(kVal-1))
        Var_qofBeta <- LofBeta^2*Var_QofBeta # Ignoring variance of within study variances
        df <- 2*((qofBeta^2)/Var_qofBeta)
  
        varDL_thetaHat <- 1/sumTrueW(kVal, tauHat_value, sigmas) # DL variance estimate
        varHK_thetaHat <- (1/((kVal-1)))*Q_RE(log_odds_ratios, thetaHat_val, kVal, betas) #HKSJ variance estimate
        
        hc1 <- hc1(log_odds_ratios, thetaHat_val, kVal, sigmas, tauHat_value) 
        hc2 <- hc2(log_odds_ratios, thetaHat_val, kVal, sigmas, tauHat_value)
        hc3 <- hc3(log_odds_ratios, thetaHat_val, kVal, sigmas, tauHat_value)
        
        #----------------------------------------------------------------------------#
        # Compute degrees of freedom adjustment for HC2 
        var_hc2_term1 <- var_hc2_term1(c_hc2, betas, sigmas, tauHat_value)
        var_hc2_term2 <- var_hc2_term2_1(c_hc2)*var_hc2_term2_2(c_hc2, betas, sigmas, tauHat_value)
        var_hc2_term3 <- var_hc2_term3_1(c_hc2, betas, sigmas, tauHat_value)*
          var_hc2_term3_2(betas, sigmas, tauHat_value)
        var_hc2_term4<- (var_hc2_term4_1_1(c_hc2)-var_hc2_term4_1_2(c_hc2))*
          var_hc2_term4_2(betas, sigmas, tauHat_value)
        var_hc2_term5 <- var_hc2_term5_1(c_hc2, betas, sigmas, tauHat_value)*
          var_hc2_term5_2(c_hc2, sigmas)
        var_hc2_term6 <- var_hc2_term6(c_hc2, betas, sigmas, tauHat_value)
        var_hc2_term7 <- var_hc2_term7(c_hc2, betas, sigmas, tauHat_value)
        var_hc2_term8 <- var_hc2_term8_1(betas, sigmas, tauHat_value)*((
          var_hc2_term8_2_1(c_hc2, betas, sigmas, tauHat_value)*var_hc2_term8_2_2(c_hc2))-
            var_hc2_term8_3(c_hc2, betas, sigmas, tauHat_value))
        
        # Compute variance of HC2
        var_hc2 <- var_hc2_term1 + var_hc2_term2 + var_hc2_term3 + var_hc2_term4 +
          var_hc2_term5 + var_hc2_term6 - var_hc2_term7 - var_hc2_term8
        
        # Compute expectetation of HC2
        mean_hc2 <- mean_hc2_term1(c_hc2, betas, sigmas, tauHat_value)+
          mean_hc2_term2_1(c_hc2, kVal)*mean_hc2_term2_2(betas, sigmas, tauHat_value) # Exp value of HC2
        
        # Compute degrees of freedom of HC2
        df_hc2_est <- 2*(mean_hc2^2)/var_hc2
        
        #-----------------------------------------------------------------------------------#
        # Compute degrees of freedom adjustment for HC3
        var_hc3_term1 <- var_hc3_term1(c_hc3, betas, sigmas, tauHat_value)
        var_hc3_term2 <- var_hc3_term2_1(c_hc3)*var_hc3_term2_2(c_hc3, betas, sigmas, tauHat_value)
        var_hc3_term3 <- var_hc3_term3_1(c_hc3, betas, sigmas, tauHat_value)*
          var_hc3_term3_2(betas, sigmas, tauHat_value)
        var_hc3_term4<- (var_hc3_term4_1_1(c_hc3)-var_hc3_term4_1_2(c_hc3))*
          var_hc3_term4_2(betas, sigmas, tauHat_value)
        var_hc3_term5 <- var_hc3_term5_1(c_hc3, betas, sigmas, tauHat_value)*
          var_hc3_term5_2(c_hc3, sigmas)
        var_hc3_term6 <- var_hc3_term6(c_hc3, betas, sigmas, tauHat_value)
        var_hc3_term7 <- var_hc3_term7(c_hc3, betas, sigmas, tauHat_value)
        var_hc3_term8 <- var_hc3_term8_1(betas, sigmas, tauHat_value)*((
          var_hc3_term8_2_1(c_hc3, betas, sigmas, tauHat_value)*var_hc3_term8_2_2(c_hc3))-
            var_hc3_term8_3(c_hc3, betas, sigmas, tauHat_value))
        
        # Compute variance of HC3
        var_hc3 <- var_hc3_term1 + var_hc3_term2 + var_hc3_term3 + var_hc3_term4 +
          var_hc3_term5 + var_hc3_term6 - var_hc3_term7 - var_hc3_term8
        
        # Compute expectetation of HC3
        mean_hc3 <- mean_hc3_term1(c_hc3, betas, sigmas, tauHat_value)+
          mean_hc3_term2_1(c_hc3, kVal)*mean_hc3_term2_2(betas, sigmas, tauHat_value) # Exp value of HC3
        
        # Compute degrees of freedom of HC3
        df_hc3_est <- 2*(mean_hc3^2)/var_hc3
        
        #-----------------------------------------------------------------------------------#
        # Compute confidence intervals for all methods
        LB_DL <- thetaHat_val-qnorm(.975)*sqrt(varDL_thetaHat)
        UB_DL <- thetaHat_val+ qnorm(.975)*sqrt(varDL_thetaHat)
        is_in_CI_DL <- ifelse((theta >= LB_DL & theta <= UB_DL), 1, 0)
        DL_CI_length <- UB_DL-LB_DL
        
        LB_HK <- thetaHat_val- qt(.975, kVal - 1)*sqrt(varHK_thetaHat) 
        UB_HK <- thetaHat_val+ qt(.975, kVal - 1)*sqrt(varHK_thetaHat) 
        is_in_CI_HK <- ifelse((theta >= LB_HK & theta <= UB_HK), 1, 0)
        HK_CI_length <- UB_HK-LB_HK
        
        LB_Hartung <- thetaHat_val- qt(.975, df)*sqrt(qofBeta) 
        UB_Hartung <- thetaHat_val+ qt(.975, df)*sqrt(qofBeta) 
        is_in_CI_Hartung <- ifelse((theta >= LB_Hartung & theta <= UB_Hartung), 1, 0)
        Hartung_CI_length <- UB_Hartung-LB_Hartung
        
        LB_hc1 <- thetaHat_val - qt(.975, kVal - 1) * sqrt(hc1)
        UB_hc1 <- thetaHat_val + qt(.975, kVal - 1) * sqrt(hc1)
        is_in_CI_hc1 <- ifelse((theta >= LB_hc1 & theta <= UB_hc1), 1, 0)
        hc1_CI_length <- UB_hc1-LB_hc1
        
        LB_hc2 <- thetaHat_val - qt(.975, kVal - 1) * sqrt(hc2)
        UB_hc2 <- thetaHat_val + qt(.975, kVal - 1) * sqrt(hc2)
        is_in_CI_hc2 <- ifelse((theta >= LB_hc2 & theta <= UB_hc2), 1, 0)
        hc2_CI_length <- UB_hc2-LB_hc2
        
        LB_hc3 <- thetaHat_val - qt(.975, kVal - 1) * sqrt(hc3)
        UB_hc3 <- thetaHat_val + qt(.975, kVal - 1) * sqrt(hc3)
        is_in_CI_hc3 <- ifelse((theta >= LB_hc3 & theta <= UB_hc3), 1, 0)
        hc3_CI_length <- UB_hc3-LB_hc3
        
        LB_hc2_newDOF_est <- thetaHat_val - qt(.975, df_hc2_est) * sqrt(hc2)
        UB_hc2_newDOF_est <- thetaHat_val + qt(.975, df_hc2_est) * sqrt(hc2)
        is_in_CI_hc2_newDOF_est <- ifelse((theta >= LB_hc2_newDOF_est & theta <= UB_hc2_newDOF_est), 1, 0)
        hc2_newDOF_est_CI_length <- UB_hc2_newDOF_est-LB_hc2_newDOF_est
        
        LB_hc3_newDOF_est <- thetaHat_val - qt(.975, df_hc3_est) * sqrt(hc3)
        UB_hc3_newDOF_est <- thetaHat_val + qt(.975, df_hc3_est) * sqrt(hc3)
        is_in_CI_hc3_newDOF_est <- ifelse((theta >= LB_hc3_newDOF_est & theta <= UB_hc3_newDOF_est), 1, 0)
        hc3_newDOF_est_CI_length <- UB_hc3_newDOF_est-LB_hc3_newDOF_est
        
        mainResults[currentRow,] <- c(kVal, currentTau, j, 
                                      df, df_hc2_est, df_hc3_est, pi_c[i], pi_t[i],
                                      thetaHat_val, tauHat_value,
                                      varDL_thetaHat, varHK_thetaHat, qofBeta, hc1, hc2, hc3,
                                      is_in_CI_DL, is_in_CI_HK, is_in_CI_Hartung,
                                      is_in_CI_hc1, is_in_CI_hc2, is_in_CI_hc3,
                                      is_in_CI_hc2_newDOF_est, is_in_CI_hc3_newDOF_est,
                                      DL_CI_length, HK_CI_length, Hartung_CI_length,
                                      hc1_CI_length, hc2_CI_length, hc3_CI_length,
                                      hc2_newDOF_est_CI_length, hc3_newDOF_est_CI_length)
        currentRow <- currentRow + 1
      }
  
    }
  }
  return(mainResults)
}

df_sims <- data.frame()
for (l in 1:10000){
  currentResults <- main()
  df_sims <- rbind(df_sims, currentResults)
}

#saveRDS(df_sims, "log-OR/oneLarge_k10_logOR.rds")

# test <- df_sims %>% 
#   filter(V7 %in% 0.1) %>% 
#   filter(V8 %in% 0.1) %>% 
#   mutate(i_sq = rep(c(0, 1/3, 1/2), 10000)) %>% 
#   group_by(V1, V3, i_sq) %>% 
#   mutate(cov_prob_DL = sum(V17)/10000,
#          cov_prob_HK = sum(V18)/10000,
#          cov_prob_HARTUNG = sum(V19)/10000,
#          cov_prob_hc1 = sum(V20)/10000,
#          cov_prob_hc2 = sum(V21)/10000,
#          cov_prob_hc3 = sum(V22)/10000,
#          cov_prob_hc2_newDOF = sum(V23)/10000,
#          cov_prob_hc3_newDOF = sum(V24)/10000)

