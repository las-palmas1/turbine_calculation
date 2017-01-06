import numpy as np


# ----------------------- Данные для расчета цикла ----------------------------

N_e = 2.0687e6
T_comb_stag = 1463.5
T_a = 288
p_a = 0.10133e6
Q_n = 43600e3
l0 = 14.61
M_v = 0
eta_comb = 0.98
sigma_comb = 0.961
eta_comp_stag_p = 0.90
eta_comp_turb_stag_p = 0.91
eta_power_turb_stag_p = 0.91
eta_m = 0.99
eta_r = 0.985
sigma_in = 0.995
sigma_out = 0.995
g_outflow = 0.01
g_cooling = 0.18
g_return = 0.02
lam_out = 0.05
pi_comp_stag = np.array([6, 9, 12, 15, 18, 21, 24, 27])
