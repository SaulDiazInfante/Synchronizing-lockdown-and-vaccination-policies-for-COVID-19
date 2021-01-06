# Title :      code based in Chatzilena2019
# Objective : TODO
# Created by: saul
# Created on: 7/1/20

library(deSolve)
library(dplyr)
library(rstan)
library(outbreaks)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#### load data from outbreaks package ####
onset <- influenza_england_1978_school$date
cases <- influenza_england_1978_school$in_bed  #Number of students in bed
write.csv(influenza_england_1978_school, "influenza_england_1978_school.csv")
N <- length(onset) # Number of days observed throughout the outbreak
pop <- 763         # Population
sample_time <- 1:N
# Modify data into a form suitable for Stan
flu_data <- list(n_obs = N,
                 n_theta = 2,
                 n_difeq = 3,
                 n_pop = pop,
                 y = cases,
                 t0 = 0,
                 ts = sample_time)

# Specify parameters to monitor
parameters <- c("y_hat", "y_init", "theta",  "R_0")
#### Model 1 - Poisson model ####
# SIR rstan model
mod1_stat <- '
functions {
  real[] SIR(real t,  // time
  real[] y,           // system state {susceptible,infected,recovered}
  real[] theta,       // parameters
  real[] x_r,
  int[] x_i) {

  real dy_dt[3];

  dy_dt[1] = - theta[1] * y[1] * y[2];
  dy_dt[2] = theta[1] * y[1] * y[2] - theta[2] * y[2];
  dy_dt[3] = theta[2] * y[2];

  return dy_dt;
  }

  }
  data {
  int<lower = 1> n_obs;       // number of days observed
  int<lower = 1> n_theta;     // number of model parameters
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower = 1> n_pop;       // population
  int y[n_obs];           // data, total number of infected individuals each day
  real t0;                // initial time point (zero)
  real ts[n_obs];         // time points observed
  }

  transformed data {
  real x_r[0];
  int x_i[0];
  }

  parameters {
  real<lower = 0> theta[n_theta]; // model parameters {beta,gamma}
  real<lower = 0, upper = 1> S0;  // initial fraction of susceptible individuals
  }

  transformed parameters{
  real y_hat[n_obs, n_difeq]; // solution from the ODE solver
  real y_init[n_difeq];     // initial conditions for both fractions of S and I

  y_init[1] = S0;
  y_init[2] = 1 - S0;
  y_init[3] = 0;
  y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);

  }

  model {
  real lambda[n_obs];      //poisson parameter

  //priors
  theta[1] ~ lognormal(0,1);
  theta[2] ~ gamma(0.004,0.02);  //Assume mean infectious period = 5 days
  S0 ~ beta(0.5, 0.5);

  //likelihood
  for (i in 1:n_obs){
  lambda[i] = y_hat[i, 2] * n_pop;
  }
  y ~ poisson(lambda);
  }

  generated quantities {
  real R_0;      // Basic reproduction number
  R_0 = theta[1]/theta[2];
  }
'
m1 <- stan_model(model_code = mod1_stat)
#### mcmcm parameters ####
n_chains <- 5
n_warmups <- 500
n_iter <- 100500
n_thin <- 50
set.seed(1234)
# Set initial values:
ini_1 <- function(){
  list(theta = c(runif(1, 0, 5), runif(1,0.2, 0.4)),
       S0 = runif(1, (pop - 3) / pop, (pop - 1) / pop))
}

time.start_nuts1 <- Sys.time()
nuts_fit_1 <- sampling(m1,
                      data = flu_data,
                      pars = parameters,
                      init = ini_1,
                      chains = n_chains,
                      warmup = n_warmups,
                      iter = n_iter,
                      thin = n_thin,
                      seed = 13219)
time.end_nuts1 <- Sys.time()
duration_nuts1 <- time.end_nuts1 - time.start_nuts1
nuts_fit_1_summary <- summary(nuts_fit_1,
                              pars = c("lp__",
                                       "theta[1]",
                                       "theta[2]",
                                       "y_init[1]",
                                       "R_0"))$summary
print(nuts_fit_1_summary, scientific = FALSE, digits = 2)
#### Post analysis ####
posts_1 <-  rstan::extract(nuts_fit_1)
mod1_diagnostics  <- rstan::get_sampler_params(nuts_fit_1)
#
# Check for divergent transitions
#
rstan::check_divergences(nuts_fit_1)
posterior_1 <- as.array(nuts_fit_1)
color_scheme_set("viridis")
#
# Markov chain traceplots
#
mcmc_trace(posterior_1, pars = "lp__")
mcmc_trace(posterior_1, pars = c("theta[1]", "theta[2]", "y_init[1]"))
mcmc_trace(posterior_1, pars = "R_0")

# Univariate and bivariate marginal posterior distributions
pairs(nuts_fit_1,
      pars = c("theta[1]", "theta[2]", "y_init[1]"),
      labels = c("beta", "gamma", "s(0)"),
      cex.labels = 1.5,
      font.labels = 9,
      condition = "accept_stat__")

# Kernel density estimates of each Markov chain separately, overlaid
mcmc_dens_overlay(posterior_1,
                  pars = c("theta[1]", "theta[2]", "y_init[1]"))

#Central posterior uncertainty intervals
mcmc_intervals(posterior_1,pars = c("theta[1]", "theta[2]", "y_init[1]"))

#### Model Fit ####
# Model fitted values across the observed time period
fit_I_1 <- posts_1$y_hat[, , 2]    # Fitted fraction of infected 
fit_SIR_1 <- fit_I_1 * pop         # Fitted number of infected
median_I_1 = apply(fit_SIR_1, 2, median)
low_I_1 = apply(fit_SIR_1, 2, quantile, probs = c(0.025))
high_I_1 = apply(fit_SIR_1, 2, quantile, probs = c(0.975))
df_sample_N = data.frame(cases, sample_time)
df_fit_I_1 = data.frame(median_I_1, low_I_1, high_I_1, sample_time)

save(df_sample_N, file = "data.Rda")
save(df_fit_I_1, file = "df_I_det_Poiss.Rda")

ggplot(df_sample_N,
       aes(x = sample_time, y = cases)) +
  geom_ribbon(aes(x = sample_time, ymin = low_I_1, ymax = high_I_1),
              fill = "orange", alpha = 0.6) +
  geom_line(data = df_fit_I_1,
            aes(x = sample_time, y = median_I_1, color = "Median"),
            size = 1.3) +
  geom_point(shape = 19, size = 3, (aes(color = "Data"))) +
  scale_colour_manual(name = '',
                      values = c('Data' = 'black',
                                 'Median' = 'darkorange3')) +
  guides(colour = guide_legend(override.aes = list(shape = c(16,NA),
                                                   linetype = c(0,1)))) +
  labs(x = "Time (days)",
       y = "Number of Infected students") +
  scale_x_continuous(limits = c(0, 14), breaks = c(0,7,14)) +
  scale_y_continuous(limits = c(0,400), breaks = c(0,100,200,300,400)) +
  theme_bw() + theme(text = element_text(size = 20))

