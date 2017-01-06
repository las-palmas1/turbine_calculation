import matplotlib.pyplot as plt
from strength_calculation.two_calculations import TwoCalculations
from strength_calculation.initial_data import *
from strength_calculation.radial_partition import RadialPartition
import strength_calculation.config as config
import os

# ---------------------------------------------
# задание исходных данных для расчета
# ---------------------------------------------
h_arr = config.h_arr
F_blade_arr = config.F_blade_arr
x_disk_end_front_arr = config.x_disk_end_front_arr
x_disk_end_back_arr = config.x_disk_end_back_arr
r_arr = config.r_arr
n = config.n
z_blade = config.z_blade
T1 = config.T1
T_m = config.T_m
material = EI698()
part_number = 20

p_m = z_blade * config.force_blade / (2 * np.pi * r_arr[len(r_arr) - 1] * h_arr[len(h_arr) - 1]) + config.pressure_tail
print(config.pressure_tail)
# ----------------------------------
# расчет на прочность
# ----------------------------------
two_calculations = TwoCalculations(n, z_blade, r_arr, h_arr, x_disk_end_front_arr, x_disk_end_back_arr, F_blade_arr,
                                   material, T1, T_m, p_m, part_number, config.d_hole_arr, config.b_hole_arr,
                                   config.r_hole_arr)


# ----------------------------------------------------
# получение массивов для построения графиков
# ----------------------------------------------------

partition = RadialPartition(n, z_blade, r_arr, h_arr, x_disk_end_front_arr, x_disk_end_back_arr, F_blade_arr,
                            material, T1, T_m, p_m, part_number)

r_arr_plot = []
h_arr_plot = []
x_disk_end_front_arr_plot = []
x_disk_end_back_arr_plot = []

for i in range(part_number + 1):
    if i == 0 or i == part_number:
        r_arr_plot.append(partition.r_arr[i])
    else:
        r_arr_plot.append(partition.r_arr[i])
        r_arr_plot.append(partition.r_arr[i])

for i in range(part_number):
    h_arr_plot.append(partition.h_arr[i])
    h_arr_plot.append(partition.h_arr[i])
    x_disk_end_front_arr_plot.append(partition.x_disk_end_front_arr[i])
    x_disk_end_front_arr_plot.append(partition.x_disk_end_front_arr[i])
    x_disk_end_back_arr_plot.append(partition.x_disk_end_back_arr[i])
    x_disk_end_back_arr_plot.append(partition.x_disk_end_back_arr[i])

sigma_r_plot = []
sigma_t_plot = []

for i in range(part_number):
    sigma_r_plot.append(two_calculations.sigma_r_arr_real[i][0])
    sigma_r_plot.append(two_calculations.sigma_r_arr_real[i][1])
    sigma_t_plot.append(two_calculations.sigma_t_arr_real[i][0])
    sigma_t_plot.append(two_calculations.sigma_t_arr_real[i][1])

# -----------------------------------------
# построение графика напряжений
# ----------------------------------------
plt.figure(figsize=(8, 6))
plt.axes([0.12, 0.15, 0.86, 0.84])
plt.xlabel(r'$\sigma,\ МПа$', fontsize=18)
plt.ylabel(r'$r,\ м$', fontsize=18)
plt.plot([item / 1e6 for item in sigma_t_plot], r_arr_plot, 'k.')
plt.plot([item / 1e6 for item in sigma_r_plot], r_arr_plot, 'k.')
plt.plot([item / 1e6 for item in two_calculations.sigma_t_arr_av], partition.r_arr, color='blue', label=r'$\sigma_t$',
         linewidth=2, linestyle='-')
plt.plot([item / 1e6 for item in two_calculations.sigma_r_arr_av], partition.r_arr, color='red', label=r'$\sigma_r$',
         linewidth=2, linestyle='-')
plt.plot([item / 1e6 for item in two_calculations.sigma_eq_arr], partition.r_arr, color='black',
         label=r'$\sigma_{экв}$', linewidth=2, linestyle='--')

# ------------------------------------------
# построение графика диска
# ------------------------------------------
dx1 = 400
sigma_r_max = max(two_calculations.sigma_r_arr_av)
sigma_t_max = max(two_calculations.sigma_t_arr_av)
x_min1 = -400
x_max1 = round(max(sigma_r_max, sigma_t_max, two_calculations.sigma_eq_max) * 1e-6, -2) + 100
k1 = 0.2 * (x_max1 - x_min1) / max(partition.h_arr)
plt.plot([item * k1 + dx1 + x_min1 for item in x_disk_end_back_arr_plot], r_arr_plot, color='black', linewidth=1)
plt.plot([item * k1 + dx1 + x_min1 for item in x_disk_end_front_arr_plot], r_arr_plot, color='black', linewidth=1)
plt.plot([x_disk_end_back_arr_plot[0] * k1 + dx1 + x_min1, x_disk_end_front_arr_plot[0] * k1 + dx1 + x_min1],
         [r_arr_plot[0], r_arr_plot[0]], color='black', linewidth=1)
plt.plot([x_disk_end_back_arr_plot[2 * part_number - 1] * k1 + dx1 + x_min1,
          x_disk_end_front_arr_plot[2 * part_number - 1] * k1 + dx1 + x_min1],
         [r_arr_plot[2 * part_number - 1], r_arr_plot[2 * part_number - 1]], color='black', linewidth=1)
plt.grid()
plt.xticks(np.linspace(x_min1, x_max1, 11))
plt.axis([x_min1, x_max1, 0,
          partition.r_arr[len(partition.r_arr) - 1] + 0.01])
plt.legend(fontsize=18, loc=0)

plt.savefig(os.path.join(os.path.dirname(__file__), 'output', 'StressPlot'))


# ---------------------------------------------
# потроение графика коэффициента запаса
# ---------------------------------------------
plt.figure(figsize=(8, 6))
plt.axes([0.12, 0.15, 0.86, 0.84])
plt.xlabel(r'$n_т$', fontsize=18)
plt.ylabel(r'$r,\ м$', fontsize=18)
plt.plot([two_calculations.safety_factor(item) for item in partition.r_arr], partition.r_arr,
         color='red', label=r'$n_т$', linewidth=2)

# ------------------------------------------
# построение графика диска
# -------------------------------------------
dx3 = 1
x_min3 = 0
x_max3 = 4
k3 = 0.2 * (x_max3 - x_min3) / max(partition.h_arr)
plt.plot([item * k3 + dx3 + x_min3 for item in x_disk_end_back_arr_plot], r_arr_plot, color='black', linewidth=1)
plt.plot([item * k3 + dx3 + x_min3 for item in x_disk_end_front_arr_plot], r_arr_plot, color='black', linewidth=1)
plt.plot([x_disk_end_back_arr_plot[0] * k3 + dx3 + x_min3, x_disk_end_front_arr_plot[0] * k3 + dx3 + x_min3],
         [r_arr_plot[0], r_arr_plot[0]], color='black', linewidth=1)
plt.plot([x_disk_end_back_arr_plot[2 * part_number - 1] * k3 + dx3 + x_min3,
          x_disk_end_front_arr_plot[2 * part_number - 1] * k3 + dx3 + x_min3],
         [r_arr_plot[2 * part_number - 1], r_arr_plot[2 * part_number - 1]], color='black', linewidth=1)
plt.grid()
plt.axis([x_min3, x_max3, 0, partition.r_arr[len(partition.r_arr) - 1] + 0.01])
plt.legend(fontsize=18, loc=1)

plt.savefig(os.path.join(os.path.dirname(__file__), 'output', 'SafetyFactorPlot'))
