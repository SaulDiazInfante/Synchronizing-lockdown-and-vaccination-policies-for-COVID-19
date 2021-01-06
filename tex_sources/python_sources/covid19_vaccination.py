import numpy as np
import json
import matplotlib.pyplot as plt
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
                                data_file_name='vaccination_data_solution2020-06-16.pkl',
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
    fig, ((ax_s, ax_e, ax_cl), (ax_i_s, ax_i_a, ax_r), (ax_d, ax_v, ax_t)) = plt.subplots(nrows=3, ncols=3)
    #
    n_whole = 1.0  # without rescaling population
    t = df['time']
    ax_s = plt.subplot(331)
    ax_s.plot(t, n_whole * df['s'], label="s")
    ax_s.legend(loc=0)
    #
    ax_e = plt.subplot(332)
    ax_e.plot(t, n_whole * df['e'], label="e")
    ax_e.legend(loc=0)
    #
    cl = n_whole * (df['s'] + df['e'] + df['i_s'] +
                    df['i_a'] + df['r'] + df['d'] +
                    df['v'] + df['treat'])
    #
    ax_cl = plt.subplot(333)
    ax_cl.plot(t, cl, label="cl")
    ax_cl.legend(loc=0)
    #
    ax_i_s = plt.subplot(334)
    ax_i_s.plot(t, n_whole * df['i_s'], label="i_s")
    ax_i_s.legend(loc=0)
    #
    ax_i_a = plt.subplot(335)
    ax_i_a.plot(t, n_whole * df['i_a'], label="i_a")
    ax_i_a.legend(loc=0)
    #
    ax_r = plt.subplot(336)
    ax_r.plot(t, n_whole * df['r'], label="r")
    ax_r.legend(loc=0)
    #
    ax_d = plt.subplot(337)
    ax_d.plot(t, n_whole * df['d'], label="d")
    ax_d.legend(loc=0)
    #
    ax_v = plt.subplot(338)
    ax_v.plot(t, n_whole * df['v'], label="v")
    ax_v.legend(loc=0)
    #
    ax_t = plt.subplot(339)
    ax_t.plot(t, n_whole * df['treat'], label="treat")
    ax_t.legend(loc=0)
    #
    plt.tight_layout()
    fig.suptitle("R0: " + str(r_zero))
    plt.savefig(fig_file_name)
    plt.show()
    return


def reproductive_number(l_0, s_0, beta_s, beta_a, p, alpha_a, alpha_s):
    r_0 = p * (p * beta_s / alpha_s + (1.0 - p) * beta_a / alpha_a) * (l_0 + s_0)
    return r_0


def rhs_vaccination(t, y, beta_s, beta_a, epsilon,
                    delta_e, delta_v, p, q,
                    alpha_a, alpha_t, alpha_s,
                    mu, mu_s, mu_a,
                    lambda_v, lambda_t):
    """"
        Version designed in May 22-2020, according to the normalized version and
        lockdown process of the above Kermack Model.
        @param beta_a:
        @param beta_s:
        @param epsilon:
        @param delta_e:
        @param delta_v:
        @param p:
        @param q:
        @param alpha_s:
        @param alpha_a:
        @param alpha_t:
        @type mu:
        @param mu_s:
        @param mu_a:
        @param lambda_v:
        @param lambda_t:
        @param y:
        @type t: object
    """
    s, e, i_s, i_a, r, d, v, treat = y
    #
    n_bar = s + e + i_s + i_a + r + v + treat
    force_infection = (beta_s * i_s + beta_a * i_a) / n_bar
    rhs_s = mu * n_bar - force_infection * s - (mu + lambda_v) * s + delta_v * v
    rhs_e = force_infection * (epsilon * v + s) - (mu + delta_e) * e
    rhs_i_s = p * delta_e * e - (mu + mu_s + alpha_s + lambda_t) * i_s - (1.0 - q) * alpha_t * treat
    rhs_i_a = (1 - p) * delta_e * e - (mu + mu_a + alpha_a) * i_a
    rhs_r = alpha_s * i_s + alpha_a * i_a + q * alpha_t * treat - mu * r
    rhs_d = mu_s * i_s + mu_a * i_a
    rhs_v = lambda_v * s - epsilon * force_infection * v - (mu + delta_v) * v
    rhs_treat = lambda_t * i_s - (mu + alpha_t) * treat
    rhs = np.array([rhs_s, rhs_e, rhs_i_s, rhs_i_a, rhs_r, rhs_d, rhs_v, rhs_treat])
    return rhs
# TODO: Modify beta parameters to estimate the exit lockdown
