# Simulates CNV_burden as negative binomial and LoF_burden as gaussian
# Determines power to detect their correlation for a grid of true correlations r

################################################################################
# Libraries
################################################################################
install.packages("MASS")
install.packages("WebPower")
install.packages("data.table")
library(MASS)
library(WebPower)
library(data.table)

################################################################################
# Parameters
################################################################################

trials = 1000    # repetitions
n = 319161       # sample size (WB)
alpha = 0.05     # significance
# parameters obtained with fitdistr
mu_1 = 1.51      # CNV_burden mean
size_1 = 0.18    # CNV_burden size 
mu_2 = 13.71     # LoF_burden mean
sd_2 = 3.78      # LoF_burden sd
corr_list <- seq(0.0005, 0.01, by = 0.0005) # Gaussian correlation
set.seed(42)

################################################################################
# STEP 1: Compute power to detect CNVburden-LoFburden across grid of r values 
################################################################################

results_df <- data.frame(
  true_corr = numeric(),
  measured_corr = numeric(),
  ci_lo = numeric(),
  ci_up = numeric(),
  power = numeric()
)

for (r in corr_list){
  print(r)
  r_estimate_lst <- numeric(trials)
  ci_lo_lst <- numeric(trials)
  ci_up_lst <- numeric(trials)
  power_lst <- numeric(trials)
  
  for (i in 1:trials){
    # CNV_burden - LoF_burden covariance matrix
    Sigma <- matrix(c(1, r, r, 1), nrow=2)
    
    # correlated gaussians from multivariate normal distribution
    gaussian_samples <- mvrnorm(n=n, mu=c(0, 0), Sigma=Sigma)
    g_1 <- gaussian_samples[,1]
    g_2 <- gaussian_samples[,2]
    
    # compute their eCDF (empirical cumulative distribution function): eCDF=rank(g)/n+1
    # this gives you p=P(X<=x)
    p_1 <- rank(g_1)/(n+1)
    p_2 <- rank(g_2)/(n+1)
    
    # compute the respective quantiles (reverse transformation of CDF) 
    # this gives the values x such that P(X<=x)=p
    x_1 <- qnbinom(p_1, size = size_1, mu = mu_1)
    x_2 <- qnorm(p_2, mean = mu_2, sd = sd_2)
    test <- cor.test(x_1, x_2)
    r_estimate <- unname(test$estimate)
    ci_lo <- unname(test$conf.int[1])
    ci_up <- unname(test$conf.int[2])
    
    power <- wp.correlation(n=n, r=r_estimate, alpha=alpha)$power
    
    r_estimate_lst[i] <- r_estimate
    ci_lo_lst[i] <- ci_lo
    ci_up_lst[i] <- ci_up
    power_lst[i] <- power
    
  }
  
  tmp_df <- data.frame(
    true_corr = r,
    measured_corr = mean(r_estimate_lst),
    ci_lo = mean(ci_lo_lst),
    ci_up = mean(ci_up_lst),
    power = mean(power_lst)
  )
  
  results_df <- rbind(results_df, tmp_df)
}


fwrite(results_df, "power_LoFburden_CNVburden_corr.csv")

