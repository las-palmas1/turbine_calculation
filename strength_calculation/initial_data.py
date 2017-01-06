from scipy.interpolate import interp1d
import math
import numpy as np


class Steel37H12H8G8MFB:
    def __init__(self):
        self.rho = 7800
        self.mu = 0.33
        self.E_1 = 1.35e11
        self.E_2 = 1.1e11
        self.T_E1 = 600
        self.T_E2 = 950
        self.alpha1 = 18e-6
        self.alpha2 = 22e-6
        self.T_alpha1 = 450
        self.T_alpha2 = 850

    def alpha(self, T):
        if T > self.T_alpha1:
            return self.alpha1 + (self.alpha2 - self.alpha1) / (self.T_alpha2 - self.T_alpha1) * (T - self.T_alpha1)
        else:
            return self.alpha1

    def E(self, T):
        if T > self.T_E1:
            return self.E_1 + (self.E_2 - self.E_1) / (self.T_E2 - self.T_E1) * (T - self.T_E1)
        else:
            return self.E_1


class Titan:
    def __init__(self):
        self.rho = 4500
        self.mu = 0.32


class TitanVT8(Titan):
    def __init__(self):
        Titan.__init__(self)
        self.T_sigma02_arr = np.array([100, 370, 470, 570, 670, 770])
        self.sigma_02_arr = np.array([840, 840, 710, 690, 610, 570]) * 1e6
        self.sigma_02_interp = interp1d(self.T_sigma02_arr, self.sigma_02_arr)

    def sigma_02(self, T):
        return self.sigma_02_interp(T)

    def E(self, T):
        return 1.15e11

    def alpha(self, T):
        return 10 * 1e-6


class TitanVT3(Titan):
    def __init__(self):
        Titan.__init__(self)
        self.T_sigma02_arr = np.array([100, 370, 470, 570, 670, 770])
        self.sigma_02_arr = np.array([760, 760, 620, 530, 500, 410]) * 1e6
        self.sigma_02_interp = interp1d(self.T_sigma02_arr, self.sigma_02_arr)
        self.T_alpha_arr = [100, 200, 300, 400, 500, 600]
        self.alpha_arr = np.array([8.6, 9.8, 10.3, 10.9, 11.4, 11.4]) * 1e-6

    def sigma_02(self, T):
        return self.sigma_02_interp(T)

    def E(self, T):
        return 1.15e11

    def alpha(self, T):
        return interp1d(self.T_alpha_arr, self.alpha_arr)(T)


class Aluminium:
    def __init__(self):
        self.rho = 2770
        self.mu = 0.34


class AluminiumAK4(Aluminium):
    def sigma_02(self, T):
        return 250 * 1e6

    def alpha(self, T):
        return 22e-6

    def E(self, T):
        return 0.72e11


class EI698:
    def __init__(self):
        self.rho = 8320
        self.mu = 0.3
        self._T_E_arr = [20, 400, 500, 600, 700, 800]
        self._E_arr = [2e11, 1.82e11, 1.75e11, 1.65e11, 1.55e11, 1.43e11]
        self._T_alpha_arr = [100, 200, 300, 400, 500, 600, 700, 800, 900]
        self._alpha_arr = [11e-6, 11.4e-6, 11.7e-6, 12.1e-6, 12.4e-6, 12.7e-6, 13.4e-6, 13.9e-6, 14.7e-6]
        self._T_sigma_02_arr = [20, 400, 500, 600, 700]
        self._sigma_02_arr = [1.22e9, 1.18e9, 1.16e9, 1.12e9, 1.04e9]

    def E(self, T):
        return interp1d(self._T_E_arr, self._E_arr)(T)

    def alpha(self, T):
        return interp1d(self._T_alpha_arr, self._alpha_arr)(T)

    def sigma_02(self, T):
        return interp1d(self._T_sigma_02_arr, self._sigma_02_arr)(T)


class InitialData:
    def __init__(self, n, z_blade, r_arr, h_arr, x_disk_end_front_arr, x_disk_end_back_arr,
                 F_blade_arr, material, T1, T_m, p_m):
        self.omega = math.pi * n / 30
        self.z_blade = z_blade
        self.r_arr_init = r_arr
        self.h_arr_init = h_arr
        self.x_disk_end_front_arr_init = x_disk_end_front_arr
        self.x_disk_end_back_arr_init = x_disk_end_back_arr
        self.F_blade_arr_init = F_blade_arr
        self.material = material
        self.T1 = T1
        self.T_m = T_m
        self.x_disk_end_back = interp1d(self.r_arr_init, self.x_disk_end_back_arr_init)
        self.x_disk_end_front = interp1d(self.r_arr_init, self.x_disk_end_front_arr_init)
        self.h = interp1d(self.r_arr_init, self.h_arr_init)
        self.F_blade = interp1d(self.r_arr_init, self.F_blade_arr_init)
        self.p_m = p_m
