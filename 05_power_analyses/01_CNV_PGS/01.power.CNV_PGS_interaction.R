# Computes power to detect CNV-PGS interaction

################################################################################
# Libraries
################################################################################
library(ggplot2)

################################################################################
# Parameteres
################################################################################

m = 1000         # repetitions
n = 96716        # sample size
alpha = 0.05/119 # significance
q_lst <- c(1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2) # CNV frequency
c_lst <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1) # interaction effect size

################################################################################
# STEP 1: Simulate CNV and PGS and compute fraction of significant effect sizes
################################################################################

results_df <- data.frame(
  CNV_freq = numeric(),
  interaction_beta = numeric(),
  power = numeric())

for (q_i in q_lst){
  for (c_j in c_lst){
    print(paste(q_i, c_j))
    counter <- 0
    for (k in 1:m){
      # variables
      CNV <- rbinom(n, 1, q_i); if (sum(CNV)==0) next
      PGS <- rnorm(n)
      d <- sqrt(1-c_j^2 * var(CNV*PGS))
      e <- d*rnorm(n)

      # model pheno
      y <- c_j*CNV*PGS + e
      model <- lm(y~CNV:PGS)
      
      # get significance of interaction effect size
      p_val <- summary(model)$coefficients['CNV:PGS','Pr(>|t|)']

      # count significant interactions
      if (p_val < alpha){
        counter <- counter + 1
      } 
    }

    # power = significant/total effect sizes
    power <- counter/m
    current_df <- data.frame(CNV_freq=q_i, 
                             interaction_beta=c_j, 
                             power=power)
    
    results_df <- rbind(results_df, current_df)
  }
}

################################################################################
# STEP 2: Plot power as function of CNV_freq and interaction effect size
################################################################################

ggplot(results_df, aes(x = factor(interaction_beta),
                       y = factor(CNV_freq),
                       fill = power)) +
  geom_tile() +
  scale_fill_gradient(low = "steelblue", high = "tomato", name = "Power") +
  geom_text(aes(label = power), color = "white", size = 3.5) +
  labs(
    x = "Interaction Effect Size",
    y = "CNV Frequency",
    title = "Power Heatmap"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

