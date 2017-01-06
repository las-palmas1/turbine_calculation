import numpy as np


r_arr = [0e-3, 47.15e-3, 68.85e-3, 118.86e-3, 122.86e-3, 126.86e-3]
h_arr = [16.43e-3, 13.14e-3, 12.22e-3, 8.73e-3, 15.87e-3, 15.87e-3]
F_blade_arr = [0 for _ in r_arr]

x_disk_end_front_arr = [-0.5 * i for i in h_arr]
x_disk_end_back_arr = [0.5 * i for i in h_arr]

T1 = 200
T_m = 700

n = 18e3

force_blade = 39.96e3
z_blade = 48

pressure_tail = 46.5e6

d_hole_arr = [8.5e-3, 5e-3]
b_hole_arr = [25.5e-3, 23e-3]
r_hole_arr = [62.5e-3, 51.65e-3]