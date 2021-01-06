functions { real[]
            LSEIR(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real dy_dt[7];
    real delta_e = 0.19607843137254904;
    real mu_s = 0.125;
    real alpha_s = 0.125;
    real N_star = y[1] + y[2] + y[3] + y[4] + y[5];
    //
    dy_dt[1]   = -1.0 * (theta[1] * y[4] + theta[2] * y[3]) * (y[1] / N_star);
    dy_dt[2]   =  (theta[1] * y[4] + theta[2] * y[3]) * (y[1] / N_star)
              - delta_e * y[2];
    dy_dt[3]  =   theta[3] * delta_e * y[2] - (theta[4] + mu_s) * y[3];
    dy_dt[4]  =  (1.0 - theta[3]) * delta_e * y[2] - theta[5] * y[4];
    dy_dt[5]   =   theta[5] * y[4] + theta[4] * y[3];
    dy_dt[6]   =   mu_s * y[3];
    dy_dt[7] =   theta[3] * delta_e * y[2];
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
    real <lower = 0> theta[n_theta]; // model parameters {beta,gamma}
    real <lower = 0, upper = 1> E0;  // initial fraction of exposed individuals
    real <lower=0, upper=1> Is0;
    real <lower=0, upper=1> Ia0;
}

transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];     // initial conditions for both fractions of S and I

    y_init[1] = n_pop - (Is0 + Ia0 + E0);
    y_init[2] = E0;
    y_init[3] = Is0;
    y_init[4] = Ia0;
    y_init[5] = 0;
    y_init[6] = 0;
    y_init[7] = Is0;

    y_hat = integrate_ode_rk45(LSEIR, y_init, t0, ts, theta, x_r, x_i);
}

model {
    real lambda[n_obs];      //poisson parameter
    //priors
    theta[1] ~ normal( 1.0, 0.3);
    theta[2] ~ normal(1.0, 0.3);
    theta[3] ~ uniform(0, 0.5);
    theta[4] ~ gamma(10, 100);
    theta[5] ~ gamma(10, 50);
    E0 ~ uniform(0, 1.104167e-05);
    Is0 ~ uniform(0, 3.312501e-06);
    Ia0 ~ uniform(0, 1.104167e-05);

  //likelihood
    for (i in 1:n_obs){
       lambda[i] = y_hat[i, 7] * n_pop;
    }
    y ~ poisson(lambda);
}

// TODO: code the corret R0
// generated quantities {
//    real R_0;      // Basic reproduction number
//    R_0 = theta[1]/theta[2];
//}
