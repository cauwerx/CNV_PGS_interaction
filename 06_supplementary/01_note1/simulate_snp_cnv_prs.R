#!/usr/bin/env Rscript

## Simulate SNP/CNV data, train SNP effects in CNV-free controls, and test
## whether held-out PRS correlates with CNV carrier status. Gaussian noise is
## scaled so the target total phenotype variance is 1.

RNGkind("L'Ecuyer-CMRG")
set.seed(20260702)

n_individuals <- 450000L
n_snps <- 300L
n_cnvs <- 5L
n_repeats <- 100L
snp_h2 <- 0.10

n_cores_env <- Sys.getenv("N_CORES", unset = NA_character_)
n_cores <- if (is.na(n_cores_env)) {
  max(1L, parallel::detectCores(logical = TRUE))
} else {
  as.integer(n_cores_env)
}
n_cores <- max(1L, min(n_cores, n_repeats))

out_dir <- "outputs"
results_file <- file.path(out_dir, "cnv_prs_correlations.csv")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

simulate_repeat <- function(r) {
  message(
    "Starting repeat ", r, " / ", n_repeats,
    " on pid ", Sys.getpid()
  )

  snp_af <- rbeta(n_snps, shape1 = 0.6, shape2 = 0.6)

  G <- matrix(0L, nrow = n_individuals, ncol = n_snps)
  for (i in seq_len(n_snps)) {
    G[, i] <- rbinom(n_individuals, size = 2L, prob = snp_af[i])
  }
  colnames(G) <- paste0("SNP", seq_len(n_snps))

  ## Scale Gaussian SNP effects so the total variance explained by SNPs is 0.1.
  ## For an additive 0/1/2 SNP, Var(G_i) = 2p_i(1 - p_i).
  snp_beta_raw <- rnorm(n_snps)
  snp_var <- 2 * snp_af * (1 - snp_af)
  snp_beta <- snp_beta_raw * sqrt(snp_h2 / sum((snp_beta_raw ^ 2) * snp_var))

  cnv_freq <- runif(n_cnvs, min = 5e-5, max = 1e-3)
  cnv_alpha <- runif(n_cnvs, min = 0.5, max = 1.5)

  C <- matrix(0L, nrow = n_individuals, ncol = n_cnvs)
  for (j in seq_len(n_cnvs)) {
    C[, j] <- rbinom(n_individuals, size = 1L, prob = cnv_freq[j])
  }
  colnames(C) <- paste0("CNV", seq_len(n_cnvs))

  genetic_value <- as.numeric(G %*% snp_beta + C %*% cnv_alpha)
  genetic_var <- var(genetic_value)
  noise_var <- 1 - genetic_var

  if (noise_var <= 0) {
    stop(
      "The genetic component variance is already >= 1; cannot add noise ",
      "while keeping total phenotype variance equal to 1."
    )
  }

  ## Use Gaussian noise scaled to the residual variance needed for target Var(Y)=1.
  noise <- rnorm(n_individuals)
  noise <- noise - mean(noise)
  noise <- as.numeric(noise / sd(noise) * sqrt(noise_var))

  Y <- genetic_value + noise

  cnv_carrier <- rowSums(C) > 0L
  non_carriers <- which(!cnv_carrier)
  carriers <- which(cnv_carrier)

  n_set_a <- floor((2 / 3) * length(non_carriers))
  set_a <- sample(non_carriers, size = n_set_a, replace = FALSE)
  set_b <- c(carriers, setdiff(non_carriers, set_a))

  y_a <- Y[set_a]
  beta_hat <- numeric(n_snps)

  for (i in seq_len(n_snps)) {
    ## Requirement: estimate each SNP effect with a single univariable lm().
    beta_hat[i] <- coef(lm(y_a ~ G[set_a, i]))[2L]
  }

  valid_snps <- which(is.finite(beta_hat))
  prs_b <- as.numeric(G[set_b, valid_snps, drop = FALSE] %*% beta_hat[valid_snps])

  repeat_cor <- rep(NA_real_, n_cnvs)
  names(repeat_cor) <- paste0("CNV", seq_len(n_cnvs))

  for (j in seq_len(n_cnvs)) {
    cnv_b <- C[set_b, j]
    repeat_cor[j] <- if (sd(cnv_b) == 0) {
      NA_real_
    } else {
      cor(prs_b, cnv_b)
    }
  }

  message(
    "Finished repeat ", r, " / ", n_repeats,
    " on pid ", Sys.getpid()
  )

  repeat_cor
}

message(
  "Running ", n_repeats, " repeats using ", n_cores, " worker process(es). ",
  "Set N_CORES to override this."
)
if (.Platform$OS.type == "windows" && n_cores > 1L) {
  warning("parallel::mclapply uses one core on Windows; running serially.")
  n_cores <- 1L
}

cor_list <- parallel::mclapply(
  seq_len(n_repeats),
  simulate_repeat,
  mc.cores = n_cores,
  mc.preschedule = FALSE,
  mc.silent = FALSE
)
cor_estimates <- do.call(rbind, cor_list)
colnames(cor_estimates) <- paste0("CNV", seq_len(n_cnvs))

cor_df <- data.frame(
  repeat_id = seq_len(n_repeats),
  cor_estimates,
  check.names = FALSE
)
write.csv(cor_df, results_file, row.names = FALSE)

message("Done.")
message("Saved correlations to: ", normalizePath(results_file))
