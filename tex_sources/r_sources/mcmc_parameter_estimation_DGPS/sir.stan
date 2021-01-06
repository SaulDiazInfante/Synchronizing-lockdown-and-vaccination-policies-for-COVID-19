functions {
  real[] SIR(real t,  // time
             real[] y,
             real[] theta,
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
    lambda[i] = y_hat[i,2]*n_pop;
  }
  y ~ poisson(lambda);
}

generated quantities {
  real R_0;      // Basic reproduction number
  R_0 = theta[1]/theta[2];
}
