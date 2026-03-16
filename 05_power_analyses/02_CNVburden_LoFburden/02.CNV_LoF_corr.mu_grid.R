# Simulates CNV_burden as negative binomial and LoF_burden as gaussian
# Determines power to detect their correlation across grid of CNVburden mean & LoFburden mean

#############################################################################################
# Libraries
#############################################################################################
install.packages("MASS")
install.packages("WebPower")
install.packages("data.table")
library(MASS)
library(WebPower)
library(data.table)

#############################################################################################
# Parameters
#############################################################################################
trials = 100                  # repetitions
n = 319161                    # sample size  (WB)
alpha = 0.05                  # significance
mu_1_list = seq(0.1, 1, 0.1)  # CNV_burden mean
# mu_1_list = seq(1, 10, 1) # CNV_burden mean
size_1 = 0.18                 # CNV_burden size 
mu_2_list = seq(5, 25, 1)     # LoF_burden mean
sd_2 = 3.78                   # LoF_burden sd
r = 1

set.seed(42)

#############################################################################################
# STEP 1: Compute power to detect CNVburden-LoFburden across grid of mu values
#############################################################################################

results_df <- data.frame(
  mu_1 = numeric(),
  mu_2 = numeric(),
  true_corr = numeric(),
  measured_corr = numeric(),
  power = numeric()
)

for (mu_1 in mu_1_list){
  for (mu_2 in mu_2_list){
    print(paste(mu_1, mu_2))
    r_estimate_lst <- numeric(trials)
    power_lst <- numeric(trials)
    
    for (i in 1:trials){
      # CNV_burden - LoF_burden covariance matrix
      Sigma <- matrix(c(1, r, r, 1), nrow=2)
    
      # correlated gaussians from multivariate normal distribution
      gaussian_samples <- mvrnorm(n=n, mu=c(0, 0), Sigma=Sigma)
      g_1 <- gaussian_samples[,1]
      g_2 <- gaussian_samples[,2]
    
      # compute their eCDF (empirical cumulative distribution function): eCDF=rank(g)/(n+1)
      # this gives you p=P(X<=x)
      p_1 <- rank(g_1)/(n+1)
      p_2 <- rank(g_2)/(n+1)
    
      # compute the neg binomial and normal quantiles (reverse transformation of CDF) 
      # this gives the values x such that P(X<=x)=p
      x_1 <- qnbinom(p_1, size=size_1, mu=mu_1)
      x_2 <- qnorm(p_2, mean=mu_2, sd=sd_2)
      
      r_estimate <- cor(x_1, x_2)
      
      power <- wp.correlation(n=n, r=r_estimate, alpha=alpha)$power
    
      r_estimate_lst[i] <- r_estimate
      power_lst[i] <- power
      
    }
  
    tmp_df <- data.frame(
      mu_1 = mu_1,
      mu_2 = mu_2,
      true_corr = r,
      measured_corr = mean(r_estimate_lst),
      power = mean(power_lst)
    )
    
    results_df <- rbind(results_df, tmp_df)
  }
}


fwrite(results_df, paste0("power_LoFburden_CNVburden_corr.mu1_mu2_grid.csv"))

