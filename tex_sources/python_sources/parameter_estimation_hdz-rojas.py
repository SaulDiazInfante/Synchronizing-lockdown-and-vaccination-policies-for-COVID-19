import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cparrarojas_example import PDEparams as pde


def LotkaVolterra(z, t, a, b):
    """
    The input z corresponds to the current state of the system, z = [x, y]. Since the input is in 1D, no
    pre-processing is needed.
    t is the current time.
    a and b correspond to the unknown parameters.
    """
    x, y = z
    return [a * x * (1 - x) - b * x * y, b * x * y - y]


def initial_x():
    return 0.3


def initial_y():
    return 0.5


str_path = '/home/saul/sauld@cimat.mx/UNISON/Articles/COVID19-VACINATION/python_sources/cparrarojas_example/examples/'

str_path = str_path + 'LotkaVolterraData.csv'
df = pd.read_csv(str_path)
df.head()
my_model = pde.PDEmodel(df,
                        LotkaVolterra,
                        [initial_x, initial_y],
                        bounds=[(2, 4), (0.5, 2)],
                        param_names=[r'$a$', r'$b$'],
                        nvars=2,
                        ndims=0,
                        nreplicates=3,
                        obsidx=None,
                        outfunc=None)

var_ic = my_model.initial_condition
my_model.fit()
var_bp = my_model.best_params
var_be = my_model.best_error
my_model.likelihood_profiles()
var_p = my_model.result_profiles
my_model.plot_profiles()
my_model.bootstrap()
var_bs = my_model.bootstrap_summary
var_br = my_model.bootstrap_raw
my_model.plot_bootstrap()
