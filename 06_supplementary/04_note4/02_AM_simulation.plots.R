#install.packages("plot3D")
#install.packages("minpack.lm")
#install.packages("Cairo")

library(ggplot2)
library(plot3D)
library(minpack.lm)
library(Cairo)


setwd("/Users/caterinacevallos/Documents/Lausanne/thesis/penetrance_at_scale/data/55_simulations/")

beta_value <- 0

################################################################################
# 1. Load data
################################################################################

long_df <- read.csv(paste0('simulations.beta_', beta_value,'.csv'))
alpha <- subset(long_df, z_axis=="alpha_beta_gamma")$alpha
gamma <- subset(long_df, z_axis=="alpha_beta_gamma")$gamma
# get product: alpha*beta*gamma
product <- subset(long_df, z_axis=="alpha_beta_gamma")$z_values
# get correlation: r(G,C)
r <- subset(long_df, z_axis=="r_GC")$z_values
# merge in compact format
df <- data.frame(
  alpha = alpha,
  gamma = gamma,
  product = product,
  r = r
)

################################################################################
# 2. Fit non-linear relationships with linear model
################################################################################

# fit (rotated) hyperbolic paraboloid (non-linear relationship):
# for the product
fit_product_lm <- lm(
  product ~ I(alpha^2) + I(gamma^2) + alpha * gamma,
  data = df
)

# for the r(G, C)
fit_r_lm <- lm(
  r ~ I(alpha^2) + I(gamma^2) + alpha * gamma,
  data = df
)

################################################################################
# 3. Predict product and r(C, G) on new grid with non-linear models
################################################################################

# new grid
grid.lines <- 80
alpha.pred <- seq(min(alpha), max(alpha), length.out = grid.lines)
gamma.pred <- seq(min(gamma), max(gamma), length.out = grid.lines)
grid <- expand.grid(alpha = alpha.pred, gamma = gamma.pred)

# make predictions on the new grid
grid$product <- predict(fit_product_lm, newdata = grid)
grid$r <- predict(fit_r_lm, newdata = grid)

product.pred <- matrix(grid$product, nrow = grid.lines, ncol = grid.lines)
r.pred <- matrix(grid$r, nrow = grid.lines, ncol = grid.lines)

################################################################################
# 4. Plot
################################################################################

CairoPDF(paste0("simulations.beta_", beta_value, ".pdf"), family = "Arial")

# scatter plot with fitted model
scatter3D(alpha, gamma, product, 
          zlim = c(-1, 1),
          # points format
          pch = 19, col="brown2", cex = 0, # cex=0 to avoid plotting scatterplot
          # ticks format
          ticktype = "detailed", nticks=2, cex.axis = 0.8,
          # surface format
          surf = list(x = alpha.pred, y = gamma.pred, z = product.pred,  facets = NA),
          # box format, colvar (color by z axis or not)
          bty = "b2", phi = 20, colvar=NULL, 
          # labs + title format 
          xlab="α", ylab='γ', zlab='', cex.lab = 1.2,  main=paste0('β=', beta_value), cex.main = 1, font.main = 1)

scatter3D(alpha, gamma, r, 
          # points format
          pch = 19, col="royalblue3", cex = 0, ticktype = "detailed",
          # surface format
          surf = list(x = alpha.pred, y = gamma.pred, z = r.pred,  facets = NA), 
          add = TRUE)

# allow drawing the legend anywhere
par(xpd = TRUE)
# display the legend for z-axis
legend(x = -0.69, y = 0.18, legend = c('α*β*γ', "r(PGS,CNV)"),
       pch = 19, col = c("brown2", "royalblue3"),
       bty = "n", cex = 0.9)

dev.off()
