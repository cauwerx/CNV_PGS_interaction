# Computes Spearman correlation between AM path coefficients (alpha, beta, gamma)

library(readxl)

DIR <- "/penetrance_at_scale/tables/"

columns <- c('ID', 'phenotype', 'cytogenic_band', 'top_model', 'alpha_r', 'alpha_SE', 'alpha_p', 'alpha_n', 'beta_r', 'beta_SE', 'beta_p', 'beta_n', 'gamma_r', 'gamma_SE', 'gamma_p', 'gamma_n', 'indirect_r', 'indirect_SE', 'total_r', 'total_SE', 'total_p', 'total_n')

df <- read_xlsx(path=paste0(DIR, 'SupTables_v3-20260305.xlsx'), sheet="TableS11_new", skip=3, col_names=columns)

df[, 5:ncol(df)] <- lapply(df[, 5:ncol(df)], as.numeric)

cat('Number of CNV-trait pairs:', nrow(df), '\n')

a_b.test <- cor.test(df$alpha_r, df$beta_r, method='spearman', exact = FALSE)
b_g.test <- cor.test(df$beta_r,  df$gamma_r, method='spearman', exact = FALSE)
a_g.test <- cor.test(df$alpha_r, df$gamma_r, method='spearman', exact = FALSE)

cat('alpha vs  beta: ', 'r=', a_b.test$estimate, 'p=', a_b.test$p.value, '\n')
cat('beta vs  gamma: ', 'r=', b_g.test$estimate, 'p=', b_g.test$p.value, '\n')
cat('alpha vs gamma: ', 'r=', a_g.test$estimate, 'p=', a_g.test$p.value, '\n')
