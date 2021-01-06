from covid19_lockdown_vaccination import *
from scipy.integrate import solve_ivp

prm = load_parameters("vaccination_ver_DPSS.json")
T = prm["T"]
n_whole = prm["n_whole"]
s_zero = prm["s_0"]
e_zero = prm["e_0"]
i_s_zero = prm["i_s_0"]
i_a_zero = prm["i_a_0"]
m_zero = prm["m_0"]
h_zero = prm["h_0"]
r_zero = prm["r_0"]
d_zero = prm["d_0"]
v_zero = prm["v_0"]
# ---------------------------------------------------------------------------------
t_eval = np.linspace(0, T, num=10000, endpoint=True)
#
z_vaccination_0 = np.array([0.0,
                            s_zero,
                            e_zero,
                            i_s_zero,
                            i_a_zero,
                            m_zero,
                            h_zero,
                            r_zero,
                            d_zero,
                            v_zero])
#
args = (prm["beta_s"],
        prm["beta_a"],
        prm["delta_h"],
        prm["kappa"],
        prm["epsilon"],
        prm["gamma_s"],
        prm["gamma_a"],
        prm["gamma_m"],
        prm["gamma_h"],
        prm["gamma_v"],
        prm["mu_i_s"],
        prm["mu_h"],
        prm["mu"],
        prm["p"],
        prm["q"],
        prm["q_hat"]
        )
#
#
sol = solve_ivp(rhs_vaccination,
                [0.0, T],
                z_vaccination_0,
                dense_output=False,
                method='LSODA',
                t_eval=t_eval,
                events=None,
                vectorized=False,
                args=args
                )
"""
kwargs = {"s_0": prm["s_0"],
          "beta_s": prm["beta_s"],
          "beta_a": prm["beta_a"],
          "p": prm["p"],
          "alpha_a": prm["alpha_a"],
          "alpha_s": prm["alpha_s"]
          }
r_zero = reproductive_number(**kwargs)
"""
df_names = ['time', 'l', 's', 'e', 'i_s', 'i_a', 'm', 'h', 'r', 'd', 'v', ]
y = np.c_[sol.t, sol.y.T]
data_solution_save(y, df_names)
vaccination_population_plot(r_zero)
# print("r_zero: ", r_zero)
# TODO: Check parameters
