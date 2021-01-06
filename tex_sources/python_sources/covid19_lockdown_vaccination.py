import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from datetime import datetime


def load_parameters(file_name='exit_lockdown_parameters.json'):
    """ Load the parameters given a JSON file.
    Parameters
    ----------
    file_name: filename with parameters values in JSON format.

    Returns
    -------
    prm: dictionary
    """
    with open(file_name) as json_file:
        prm = json.load(json_file)
    return prm


def data_solution_save(sol, column_names, file_name_prefix='vaccination_data_solution'):
    """ Save scipy.integrate.solve_ivp output as pandas DataFrame
    Parameters
    ----------
    sol: dict_values.
    column_names: Names for the columns of the DataFrame.
    file_name_prefix: filename  of desire DataFrame file.

    Returns
    -------
    df.plk: a pkl pandas DataFrame file.
    """
    df = pd.DataFrame(sol)
    df.columns = column_names
    string_date = str(datetime.date(datetime.now()))
    file_name = file_name_prefix + string_date + ".pkl"  #
    df.to_pickle(file_name)
    return


def vaccination_population_plot(r_zero,
                                data_file_name='vaccination_data_solution2020-06-20.pkl',
                                fig_file_name='model_solution.pdf'):
    """
    Plot the output solution given the pandas data frame
    output from scipy.integrate.solve_ivp
    ----------
    file_name: filename with pandas data frame  output from scipy.integrate.solve_ivp

    Returns
    -------
    Shows figures in the default backend and return a pdf figure file
    """
    df = pd.read_pickle(data_file_name)
    fig = plt.figure(constrained_layout=True)
    spec = gridspec.GridSpec(ncols=4, nrows=3, figure=fig)
    #
    ax_l = fig.add_subplot(spec[0, 0])
    ax_s = fig.add_subplot(spec[0, 1])
    ax_e = fig.add_subplot(spec[0, 2])
    ax_i_s = fig.add_subplot(spec[0, 3])
    #
    ax_i_a = fig.add_subplot(spec[1, 0])
    ax_m = fig.add_subplot(spec[1, 1])
    ax_h = fig.add_subplot(spec[1, 2])
    ax_r = fig.add_subplot(spec[1, 3])
    #
    ax_d = fig.add_subplot(spec[2, 0])
    ax_v = fig.add_subplot(spec[2, 1])
    ax_cl = fig.add_subplot(spec[2, 2:])
    #
    n_whole = 1.0  # without rescaling population
    t = df['time']
    ax_l.plot(t, n_whole * df['l'], label="l")
    ax_l.legend(loc=0)
    #
    ax_s.plot(t, n_whole * df['s'], label="s")
    ax_s.legend(loc=0)
    #
    ax_e.plot(t, n_whole * df['e'], label="e")
    ax_e.legend(loc=0)

    ax_i_s.plot(t, n_whole * df['i_s'], label="i_s")
    ax_i_s.legend(loc=0)
    #
    ax_i_a.plot(t, n_whole * df['i_a'], label="i_a")
    ax_i_a.legend(loc=0)
    #
    ax_m.plot(t, n_whole * df['m'], label="m")
    ax_m.legend(loc=0)
    #
    ax_h.plot(t, n_whole * df['h'], label="h")
    ax_h.legend(loc=0)
    #
    ax_r.plot(t, n_whole * df['r'], label="r")
    ax_r.legend(loc=0)
    #
    ax_d.plot(t, n_whole * df['d'], label="d")
    ax_d.legend(loc=0)
    #
    ax_v.plot(t, n_whole * df['v'], label="v")
    ax_v.legend(loc=0)
    #
    cl = n_whole * (df['l'] + df['s'] + df['e'] +
                    df['i_s'] + df['i_a'] + df['m'] +
                    df['h'] + df['r'] + df['d'] +
                    df['v'])
    #
    ax_cl.plot(t, cl, label="cl")
    ax_cl.legend(loc=0)
    #
    plt.tight_layout()
    fig.suptitle("R0: " + str(r_zero))
    plt.savefig(fig_file_name)
    plt.show()
    return


def reproductive_number(l_0, s_0, beta_s, beta_a, p, alpha_a, alpha_s):
    r_0 = p * (p * beta_s / alpha_s + (1.0 - p) * beta_a / alpha_a) * (l_0 + s_0)
    return r_0


def rhs_vaccination(t, y, beta_s, beta_a,
                    delta_h, kappa, epsilon,
                    gamma_s, gamma_a,
                    gamma_m, gamma_h,
                    gamma_v,
                    mu_i_s, mu_h, mu,
                    p, q, q_hat):
    """"
        Version designed in May 22-2020, according to the normalized version and
        lockdown process of the above Kermack Model.
        @param t: integration interval time
        @param y: vector of epidemiological classes
        @param beta_a: asymptomatic rate transmission
        @param beta_s: symptomatic rate transmission
        @param delta_h: transition rate to hospitalized class
        @param kappa: the inverse of mean time in exposed class
        @param epsilon: disease transmission rate in lockdown class
        @param gamma_s: the inverse of mean recover time of symptomatic class
        @param gamma_a: the inverse of mean recover time of asymptomatic class
        @param gamma_m: the inverse of mean recover time of symptomatic medicated
        @param gamma_h: the inverse of mean recover time of symptomatic hospitalized
        @param gamma_v: the inverse of mean recover time of symptomatic vaccinated
        @param mu_i_s: induced symptomatic disease mortality
        @param mu_h: per capita death rate of hospitalized class
        @param mu: per capita natural death
        @param p: fraction of exposed individuals that becomes in symptomatic case
        @param q: treatment efficiency
        @param q_hat: vaccine efficiency
    """
    l, s, e, i_s, i_a, m, h, r, d, v = y
    u_l = 0.01
    u_h = 0.01
    u_m = 0.01
    u_v = 0.000001
    # n_start = l + s + e + i_s + i_a + m + h + r + v
    n_start = s + e + i_s + i_a + m + h + r + v
    force_infection = (beta_s * i_s + beta_a * i_a) / n_start
    # rhs_l = -epsilon * force_infection * l - u_l * l - mu_l * l
    rhs_s = mu * n_start + u_l * l + (1 - q_hat) * gamma_v * v - force_infection * s - u_v * s - mu * s
    rhs_e = force_infection * (epsilon * l + s) - (kappa + mu) * e
    rhs_i_s = p * kappa * e - (gamma_s + mu_i_s + delta_h) * i_s \
              - u_m * i_s + (1 - q) * gamma_m * m - mu * i_s
    rhs_i_a = (1 - p) * kappa * e - (gamma_a + mu) * i_a
    rhs_m = u_m * i_s - (gamma_m + mu) * m
    rhs_h = delta_h * i_s - (gamma_h + mu_h) * h - (u_h + mu) * h
    rhs_r = gamma_s * i_s + gamma_a * i_a + gamma_h * h + q * gamma_m * m + q_hat * gamma_v * v + u_h * h - mu * r
    rhs_d = mu_i_s * i_s + mu_h * h
    rhs_v = u_v * s - (mu + gamma_v) * v
    rhs = np.array([l, rhs_s, rhs_e, rhs_i_s, rhs_i_a, rhs_m, rhs_h, rhs_r, rhs_d, rhs_v])
    return rhs
# TODO: Modify beta parameters to estimate the exit lockdown
