# Total vs indirect CNV-PGS correlations in mating simulations

library(tidyr)

partner_choice <- function(gamma, beta, alpha, N, n_mate, simN) {
  corz = rep(0,simN) # vector to save results
  for (j in 1:simN){
    C = rbinom(N,1,.5) # CNV ~ Bernoulli 
    tmp = var(C)*(1-alpha^2)/alpha^2 
    Y = alpha*C+sqrt(tmp)*rnorm(N) # pheno of partner 1 is function of CNV + error term
    
    G = rep(0,N)
    for (i in 1:N){
      g = rnorm(n_mate) # normal distribution of PGS = partner pool
      X = gamma*g+sqrt(1-gamma^2)*rnorm(n_mate) # phenotypes of possible mates are a function of PGS + error term
      y = Y[i] # select first partner 1 iteratively
      x = rnorm(1,mean=beta*y,sd = sqrt(1-beta^2)) # get the ideal phenotype of partner 2
      ix = which((X-x)^2==min((X-x)^2)) # from the partner pool (X) get the individual whose phenotype would have a corr=beta with phenotype of partner 1
      G[i] = g[ix[1]] # get PGS of the best partner 2 
    }
    corz[j] = cor(G,C)
  }
  return(mean(corz))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide beta value.")
}
beta <- as.numeric(args[1])

N = 1e4
n_mate = 5e2
simN = 1e2

alpha_list <- seq(-1, 1, 0.1)
gamma_list <- seq(-1, 1, 0.1)
grid <- expand.grid(alpha = alpha_list, gamma = gamma_list)

# correlation through AM
grid$alpha_beta_gamma <- grid$alpha * beta * grid$gamma
# total CNV, PGS correlation
grid$r_GC <- mapply(partner_choice, alpha=grid$alpha, gamma=grid$gamma, MoreArgs = list(beta = beta, N = N, n_mate = n_mate, simN = simN))

# add predicted values 
grid_long <- grid[,c('alpha', 'gamma', 'alpha_beta_gamma', 'r_GC')] |>
  pivot_longer(cols = c(alpha_beta_gamma, r_GC), names_to = "z_axis", values_to = "z_values")

write.csv(grid_long, paste0("simulations.beta_", beta, ".csv"), row.names = FALSE, quote = FALSE)

