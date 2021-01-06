library(deSolve)
library(dplyr)
library(rstan)
library(outbreaks)
library(bayesplot)
library(data.table)
#
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
data_path <-
    paste("/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19",
          "NovelCovid19-ControlModelling/COVID19-VACINATION/r_sources/data",
          "culiacan_prevalence_data.csv",
          sep = "/"
          )
#
covid19_culiacan_data <- fread(data_path,
                               select = c("FECHA_SINTOMAS",
                                          "i_s",
                                          "cumulative_i_s"))
covid19_culiacan_data <- data.frame(covid19_culiacan_data)
onset <- covid19_culiacan_data %>%
    select(FECHA_SINTOMAS)
#
cum_cases <- covid19_culiacan_data %>%
    select(cumulative_i_s)
cum_cases <- unlist(cum_cases, use.names = FALSE)
N <- nrow(onset)
# Number of days observed throughout the outbreak
#
pop <- 962871         # Population
sample_time <- 1:N
# Modify data into a form suitable for Stan
#
covid19_data <- list(n_obs = N,
                 n_theta = 5,
                 n_difeq = 7,
                 n_pop = pop,
                 y = cum_cases,
                 t0 = 0,
                 ts = sample_time)
# Specify parameters to monitor
parameters <- c("y_hat", "y_init", "theta")
#TODO: Recode for the correct R0
#parameters <- c("y_hat", "y_init", "theta",  "R_0")
#### Model 1 - Poisson model ####
# SIR rstan model
#
m1 <- stan_model("lockdown_mcKendrick.stan")
#### mcmcm parameters ####
n_chains <- 5
n_warmups <- 500
n_iter <- 100500
n_thin <- 50
set.seed(1234)
# Set initial values:
ini_1 <- function(){
  list(theta = c(rnorm(1, mean = 1.0, sd = 0.3),
                 rnorm(1, mean = 1.0, sd = 0.3),
                 runif(1, 0, 0.5),
                 rgamma(1, 10, 100),
                 rgamma(1, 10, 50)
                ),
       E0 = 0.0,
       Is0 = runif(1, 1, 3 / pop),
       Ia0 = runif(1, 1, 10 / pop)
       )
}
#TODO: Tune the parameter adapt_delta
#stan(..., control = list(adapt_delta = 0.99))
time.start_nuts1 <- Sys.time()
nuts_fit_1 <- sampling(m1,
                      data = covid19_data,
                      pars = parameters,
                      init = ini_1,
                      chains = n_chains,
                      warmup = n_warmups,
                      iter = n_iter,
                      thin = n_thin,
                      seed = 13219,
                      control = list(adapt_delta = 0.99))
time.end_nuts1 <- Sys.time()
duration_nuts1 <- time.end_nuts1 - time.start_nuts1
nuts_fit_1_summary <- summary(nuts_fit_1,
                              pars = c("lp__",
                                       "theta[1]",
                                       "theta[2]",
                                       "theta[3]",
                                       "theta[4]",
                                       "theta[5]",
                                       "y_init[2]",
                                       "y_init[3]",
                                       "y_init[4]"))$summary
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
mcmc_trace(posterior_1, pars = c("theta[1]",
                                 "theta[2]",
                                 "theta[3]",
                                 "theta[4]",
                                 "theta[5]",
                                 "y_init[2]",
                                 "y_init[3]",
                                 "y_init[4]"))
# mcmc_trace(posterior_1, pars = "R_0")
#
# Univariate and bivariate marginal posterior distributions
#
pairs(nuts_fit_1,
      pars = c("theta[1]",
               "theta[2]",
               "theta[3]",
               "theta[4]",
               "theta[5]",
               "y_init[2]",
               "y_init[3]",
               "y_init[4]"
               ),
      labels = c("beta_s", "beta_a", "rho"),
      cex.labels = 1.5,1
      font.labels = 9,
      condition = "accept_stat__")
#
# Kernel density estimates of each Markov chain separately, overlaid
#
mcmc_dens_overlay(posterior_1,
                  pars = c("theta[1]", "theta[2]", "y_init[1]"))
#
#Central posterior uncertainty intervals
mcmc_intervals(posterior_1 ,
               pars = c("theta[1]",
                        "theta[2]",
                        "theta[3]",
                        "theta[4]",
                        "theta[5]",
                        "y_init[2]",
                        "y_init[3]",
                        "y_init[4]",
                        ))
#
#### Model Fit ####
# Model fitted values across the observed time period
fit_I_1 <- posts_1$y_hat[, , 7]    # Fitted fraction of infected
fit_SIR_1 <- fit_I_1 * pop         # Fitted number of infected
median_I_1 = apply(fit_SIR_1, 2, median)
low_I_1 = apply(fit_SIR_1, 2, quantile, probs = c(0.025))
high_I_1 = apply(fit_SIR_1, 2, quantile, probs = c(0.975))
df_sample_N = data.frame(cum_cases, sample_time)
df_fit_I_1 = data.frame(median_I_1, low_I_1, high_I_1, sample_time)
#
save(df_sample_N, file = "data.Rda")
save(df_fit_I_1, file = "df_I_det_Poiss.Rda")
#
ggplot(df_sample_N,
       aes(x = sample_time, y = cum_cases)) +
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
    scale_x_continuous(limits = c(0, 13)) +
    scale_y_continuous(limits = c(0, 50),
                       breaks = c(0, 20, 40, 60, 80, 100)) +
  theme_bw() + theme(text = element_text(size = 20))

#### Divergence analysis
