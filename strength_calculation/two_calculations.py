from strength_calculation.radial_partition import RadialPartition
from strength_calculation.calculation import Calculation
from scipy.interpolate import interp1d
import numpy as np


class TwoCalculations(RadialPartition):
    def __init__(self, n, z_blade, r_arr, h_arr, x_disk_end_front_arr, x_disk_end_back_arr, F_blade_arr, material,
                 T1, T_m, p_m, part_number, d_hole_arr, b_hole_arr, r_hole_arr):
        RadialPartition.__init__(self, n, z_blade, r_arr, h_arr, x_disk_end_front_arr, x_disk_end_back_arr,
                                 F_blade_arr, material, T1, T_m, p_m, part_number)
        self.sigma_t11 = 100e6
        self.sigma_r11 = 100e6
        self.first_calc = Calculation(self.sigma_t11, self.sigma_r11, self.omega, self.theta, self.theta_arr, self.r_arr,
                                      self.E_arr, self.h_arr, self.rho_new_arr, self.material.mu)
        self.second_calc = Calculation(self.sigma_t11, self.sigma_r11, 0, self.theta0, self.theta0_arr, self.r_arr,
                                       self.E_arr, self.h_arr, self.rho_new_arr, self.material.mu)

        self.k = (self.p_m -
                  self.first_calc.sigma_r_arr[self.part_number -
                                              1][1]) / self.second_calc.sigma_r_arr[self.part_number - 1][1]

        self.sigma_r_arr_real = []
        self.sigma_t_arr_real = []
        self.sigma_r_arr_av = []
        self.sigma_t_arr_av = []
        self.sigma_eq_arr = []
        self.safety_factor_arr = []
        self.d_hole_arr = d_hole_arr
        self.b_hole_arr = b_hole_arr
        self.r_hole_arr = r_hole_arr

        for i in range(part_number):
            self.sigma_r_arr_real.append([])
            self.sigma_t_arr_real.append([])
            self.sigma_r_1_real = self.first_calc.sigma_r_arr[i][0] + self.k * self.second_calc.sigma_r_arr[i][0]
            self.sigma_r_2_real = self.first_calc.sigma_r_arr[i][1] + self.k * self.second_calc.sigma_r_arr[i][1]
            self.sigma_t_1_real = self.first_calc.sigma_t_arr[i][0] + self.k * self.second_calc.sigma_t_arr[i][0]
            self.sigma_t_2_real = self.first_calc.sigma_t_arr[i][1] + self.k * self.second_calc.sigma_t_arr[i][1]
            self.sigma_r_arr_real[i].append(self.sigma_r_1_real)
            self.sigma_r_arr_real[i].append(self.sigma_r_2_real)
            self.sigma_t_arr_real[i].append(self.sigma_t_1_real)
            self.sigma_t_arr_real[i].append(self.sigma_t_2_real)

        for i in range(part_number + 1):
            if i == 0:
                self.sigma_r_av = self.sigma_r_arr_real[i][0]
                self.sigma_t_av = self.sigma_t_arr_real[i][0]
            elif i != part_number:
                self.sigma_r_av = 0.5 * (self.sigma_r_arr_real[i - 1][1] + self.sigma_r_arr_real[i][0])
                self.sigma_t_av = 0.5 * (self.sigma_t_arr_real[i - 1][1] + self.sigma_t_arr_real[i][0])
            elif i == part_number:
                self.sigma_r_av = self.sigma_r_arr_real[i - 1][1]
                self.sigma_t_av = self.sigma_t_arr_real[i - 1][1]
            if self.sigma_t_av < 0:
                self.sigma_eq = np.sqrt(self.sigma_r_av ** 2 + self.sigma_t_av ** 2 - self.sigma_r_av * self.sigma_t_av)
            else:
                self.sigma_eq = self.sigma_r_av
            self.sigma_r_arr_av.append(self.sigma_r_av)
            self.sigma_t_arr_av.append(self.sigma_t_av)
            self.sigma_eq_arr.append(self.sigma_eq)
            self.safety_factor_arr.append(self.sigma_02_arr[i] / self.sigma_eq)

        self.sigma_eq_max = max(self.sigma_eq_arr)
        self.safety_factor_min = min(self.safety_factor_arr)
        self._sigma_r = interp1d(self.r_arr, self.sigma_r_arr_av)
        self._sigma_t = interp1d(self.r_arr, self.sigma_t_arr_av)
        self._sigma_eq = interp1d(self.r_arr, self.sigma_eq_arr)

        self.k_hole_arr = []
        self.sigma_t_hole = []
        self.safety_factor_hole = []
        for n, i in enumerate(self.r_hole_arr):
            k = 3 - self.d_hole_arr[n] / self.b_hole_arr[n] - self.sigma_r(i) / self.sigma_t(i)
            self.k_hole_arr.append(k)
            self.sigma_t_hole.append(k * self.sigma_t(i))
            self.safety_factor_hole.append(self.sigma_02(i) / k * self.sigma_t(i))

    def sigma_r(self, r):
        return self._sigma_r(r)

    def sigma_t(self, r):
        return self._sigma_t(r)

    def sigma_eq(self, r):
        return self._sigma_eq(r)

    def safety_factor(self, r):
        return self.sigma_02(r) / self._sigma_eq(r)
