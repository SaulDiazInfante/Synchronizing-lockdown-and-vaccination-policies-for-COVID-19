from covid19_vaccination import *
from scipy.integrate import solve_ivp

prm = load_parameters("model_vaccination_parameters.json")
T = prm["T"]
n_whole = prm["n_whole"]
s_zero = prm["s_0"]
e_zero = prm["e_0"]
i_s_zero = prm["i_s_0"]
i_a_zero = prm["i_a_0"]
r_zero = prm["r_0"]
d_zero = prm["d_0"]
v_zero = prm["v_0"]
treat_zero = prm["treat_0"]
# ---------------------------------------------------------------------------------
t_eval = np.linspace(0, T, num=10000, endpoint=True)
#
z_lockdown_0 = np.array([s_zero,
                         e_zero,
                         i_s_zero,
                         i_a_zero,
                         r_zero,
                         d_zero,
                         v_zero,
                         treat_zero
                         ])
#
args = (prm["beta_s"], prm["beta_a"], prm["epsilon"],
        prm["delta_e"], prm["delta_v"],
        prm["p"], prm["q"],
        prm["alpha_a"], prm["alpha_t"], prm["alpha_s"],
        prm["mu"], prm["mu_s"], prm["mu_a"],
        prm["lambda_v"], prm["lambda_t"])
#
#
sol = solve_ivp(rhs_vaccination,
                [0.0, T],
                z_lockdown_0,
                dense_output=False,
                method='RK45',
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
df_names = ['time', 's', 'e', 'i_s', 'i_a', 'r', 'd', 'v', 'treat']
y = np.c_[sol.t, sol.y.T]
data_solution_save(y, df_names)
vaccination_population_plot(r_zero, data_file_name='vaccination_data_solution2020-07-15.pkl')
print("r_zero: ", r_zero)
# TODO: Check parameters
