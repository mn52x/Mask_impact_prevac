rm(list = ls())

suppressPackageStartupMessages({
library(cmdstanr)
library(dplyr); library(tidyr); library(stringr)
library(ggplot2); library(posterior); library(tibble)
library(patchwork); library(here)
})

# ----- load data ----
country_name <- #country iso3 code
data <- read.csv("../country_case_mask.csv")
pop <- read.csv("../population.csv")
fix_end_df <- read.csv("../lase_observed_time_in_sim.csv")
rep_cov_df <- read.csv("../reporting_coverage.csv")

#----calculate incidence using report coverage---
#time when counterfactual starts
fix_end <- fx_end_df[fx_end_df$country == country_name, "fix_end"]  #The last date of the observed value in the simulation

#set reporting coverage 
n <- nrow(data)
idx <- 1:n
rep_cov <-  rep_cov_df[report_cov_df$country == country_name, "reporting_coverage"]
reportcoverage <- case_when(
  idx <= n ~ rep_cov 
)
#save as data column 
data <- alldata[alldata$country == country_name, ]
X0019 <- as.numeric(data$case_age_00_19)
X2039 <- as.numeric(data$case_age_20_39)
X4059 <- as.numeric(data$case_age_40_59)
X60over <- as.numeric(data$case_age_60over)

#raw observed vector
obs0019 <- pmax(data$case_age_0_19, 1e-3)
obs2039 <- pmax(data$case_age_20_39, 1e-3)
obs4059 <- pmax(data$case_age_40_59, 1e-3)
obs60 <- pmax(data$case_age_60over, 1e-3)

#incidence vector (observed/report coverage)
data$X0019_infec <- X0019 / reportcoverage
data$X2039_infec <- X2039 / reportcoverage
data$X4059_infec <- X4059 / reportcoverage
data$X60over_infec <- X60over / reportcoverage
data$total_infec <-data$X0019_infec + data$X2039_infec + data$X4059_infec + data$X60over_infec

vec0019 <- pmax(data$X0019_infec, 1e-3)
vec2039 <- pmax(data$X2039_infec, 1e-3)
vec4059 <- pmax(data$X4059_infec, 1e-3)
vec60 <- pmax(data$X60over_infec, 1e-3)

# incidence vector
vec_mat <- cbind(vec0019, vec2039, vec4059, vec60)
Y <- vec_mat

#---population : N_a-----
maching_pop <- pop[pop$country == country_name, ]
n_0019 <- maching_pop[, "case_age_0_19"]
n_2039 <- maching_pop[, "case_age_20_39"]
n_4059 <- maching_pop[, "case_age_40_59"]
n_60 <- maching_pop[, "case_age_60over"]

#---Serial interval (Nishiura, et al, 2020, weibull)----
gi_fit = list(shape=2.203339, scale=5.419875)
generation <- function(t){
  pweibull(t, shape = gi_fit$shape, scale = gi_fit$scale) -
    pweibull(t-1, shape = gi_fit$shape, scale = gi_fit$scale)
}
#-----14days discrete serial interval vector and normalize g_vec ------
g_vec <- pweibull(1:14, shape = gi_fit$shape, scale = gi_fit$scale) - 
  pweibull(0:13, shape = gi_fit$shape, scale = gi_fit$scale)
g_vec_sum <- sum(g_vec) 
g_vec_norm <- g_vec / g_vec_sum

# ------susceptible proportion------
s_mat <- matrix(NA_real_, nrow = T, ncol = 4)  

for (t in 1:T) {
  idx <- seq_len(t - 1)          
  cum_i_0019 <- if (length(idx)) sum(vec0019[idx]) else 0
  cum_i_2039 <- if (length(idx)) sum(vec2039[idx]) else 0
  cum_i_4059 <- if (length(idx)) sum(vec4059[idx]) else 0
  cum_i_60 <- if (length(idx)) sum(vec60[idx]) else 0

  s_mat[t, ] <- 1 - c(cum_i_0019/n_0019, cum_i_2039/n_2039,
                      cum_i_4059/n_4059, cum_i_60/n_60)
}

#------next generation matrix-----
#original next generation matrix from Sung-mok el al.,
A <- matrix(
  c(0.69, 0.11, 0.14, 0.04, 
    0.37, 0.87, 0.52, 0.24, 
    0.34, 0.38, 0.51, 0.22, 
    0.11, 0.22, 0.28, 0.50), 
  nrow = 4, ncol = 4, byrow = TRUE
)
eig_vals <- eigen(A %*% t(A))$values  
max_eig <- max(Re(eig_vals))
A_norm <- A/max_eig　#normalize

#step for every 28 days
#  data$date is Date 
dates <- as.Date(data$date) 
stopifnot(!anyNA(dates))
d0 <- min(dates, na.rm = TRUE)

# 0,1,2,...day, devided by 28 -> step  +1 → 1,2,3,... 
step_id <- as.integer(floor(as.numeric(dates - d0) / 28)) + 1L
n_step  <- max(step_id)
stopifnot(length(step_id) == T,
          all(step_id >= 1L),
          identical(sort(unique(step_id)), seq_len(n_step)))

#---mask-effect matrix------
maskeffect <- 0.287682072451781 #Leech et al.,
lambda_1 <- 0.02326863 #Olillia et al;, PLOS One (2022)
lambda_2 <- 0.26441345

mask0019 <- as.numeric(data$mask_age_00_19)
mask2039 <- as.numeric(data$mask_age_20_39)
mask4059 <- as.numeric(data$mask_age_40_59)
mask60 <- as.numeric(data$mask_age_60over)

#  exp(-λ * ·) for each elements
mask_w <- cbind(mask_00_19, mask_20_39, mask_40_59, mask_60over)   # T×4

mask_sus <- exp(-lambda_1 * mask_w)  # T×4
mask_inf <- exp(-lambda_2 * mask_w)  # T×4

n_age <- 4

mask_sus_mat <- array(0.0, dim = c(T, n_age, n_age))
mask_inf_mat <- array(0.0, dim = c(T, n_age, n_age))

for (t in seq_len(T)) {
  mask_sus_mat[t, , ] <- diag(mask_sus[t, ])
  mask_inf_mat[t, , ] <- diag(mask_inf[t, ])
}

#-----counterfactual-------
#scenario 1 : higher coverage
mask_CF1_0019 <- rep(0.97, T)
mask_CF1_2039 <- rep(0.97, T)
mask_CF1_4059 <- rep(0.97, T)
mask_CF1_60 <- rep(0.97, T)
mask_CF1 <- cbind(mask_CF1_0019, mask_CF1_2039, mask_CF1_4059, mask_CF1_60)  
# mask reduction for infectee
mask_sus_CF1  <- exp(-lambda_1 * mask_CF1 )  # T×4
mask_inf_CF1  <- exp(-lambda_2 * mask_CF1 ) 
# 4×4×T 
mask_sus_mat_CF1 <- array(0.0, dim = c(T, n_age, n_age))
mask_inf_mat_CF1 <- array(0.0, dim = c(T, n_age, n_age))

for (t in seq_len(T)) {
  mask_sus_mat_CF1[t, , ] <- diag(mask_sus_CF1[t, ])
  mask_inf_mat_CF1[t, , ] <- diag(mask_inf_CF1[t, ])
}

#senario 2 :lower coverage
##Asia : 1.8%
#Europe : 0.2%

mask_CF2_0019 <- rep(0.018, T)
mask_CF2_2039 <- rep(0.018, T)
mask_CF2_4059 <- rep(0.018, T)
mask_CF2_60 <- rep(0.018, T)
mask_CF2 <- cbind(mask_CF2_0019, mask_CF2_2039, mask_CF2_4059, mask_CF2_60)  
# mask reduction for infectee
mask_sus_CF2  <- exp(-lambda_1 * mask_CF2)  # T×4
mask_inf_CF2  <- exp(-lambda_2 * mask_CF2 ) 
# 4×4×T setting
mask_sus_mat_CF2 <- array(0.0, dim = c(T, n_age, n_age))
mask_inf_mat_CF2 <- array(0.0, dim = c(T, n_age, n_age))

for (t in seq_len(T)) {
  mask_sus_mat_CF2[t, , ] <- diag(mask_sus_CF2[t, ])
  mask_inf_mat_CF2[t, , ] <- diag(mask_inf_CF2[t, ])
}

#-----MCMC-------
# vector check
stopifnot(is.matrix(Y), nrow(Y) == T, ncol(Y) == 4)
stopifnot(length(step_id) == T)
stopifnot(max(step_id) == n_step)           # = step_id  max
stopifnot(all(step_id >= 1L), all(step_id <= n_step))

# Observed:
M_mat <- as.matrix(data[, c("X19", "X2039", "X4059", "X60over")])
storage.mode(M_mat) <- "integer" 
T <- nrow(data)
n_age <- 4
  
# cumlative period(ex: t=80〜T）
cum_start <- 1L
cum_end   <- T                             
if (cum_end < cum_start) { tmp <- cum_start; cum_start <- cum_end; cum_end <- tmp }
stopifnot(cum_start >= 1L, cum_end <= T)

# time flag for liklihood
use_lik_t_vec <- rep(1L, T)

#data for Stan
stan_data <- list(
  T = T,
  n_age = 4,
  n_step = n_step,                  # = max(step_id)
  step_id = as.integer(step_id),
 
  M = M_mat,  # Stan: array[T, n_age] int M
  
  # Stan: matrix[T, n_age] Y / vec_mat / s_mat
  Y = Y,
  vec_mat = vec_mat,
  s_mat = s_mat,
  
  A = A_norm,
  s_mat_init = s_mat[1:2, ],        # inisial values for t = 1,2
  
  reportcoverage = as.vector(reportcoverage),  # vector[T] (0~1)
  
  g_len = length(g_vec_norm),
  g_vec = as.vector(g_vec_norm),
  
  #population 
  n_pop = c(n_0019, n_2039, n_4059, n_60),
  
  #mask-effect diagonal matrix
  mask_sus_mat = mask_sus_mat,
  mask_inf_mat = mask_inf_mat,
  mask_sus_mat_CF1 = mask_sus_mat_CF1,
  mask_inf_mat_CF1 = mask_inf_mat_CF1,
  mask_sus_mat_CF2 = mask_sus_mat_CF2,
  mask_inf_mat_CF2 = mask_inf_mat_CF2,
  
  #start time of counterfactual
  fix_end = fix_end,
  
  # cumlative 
  cum_start = cum_start,
  cum_end   = cum_end,
)


str(stan_data)

# .stan
mod <- cmdstan_model("maskage_Poisson_log.stan")
# sampling setting
fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,              
  parallel_chains = 4,
  iter_warmup   = 500,　#for test    
  iter_sampling = 1000, # for test
  adapt_delta   = 0.96,     
  max_treedepth = 15,       
  refresh = 200,
  save_cmdstan_config = TRUE
)

#save
fit$save_output_files()
