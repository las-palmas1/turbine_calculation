from strength_calculation.initial_data import InitialData
import numpy as np


class RadialPartition(InitialData):
    def __init__(self, n, z_blade, r_arr, h_arr, x_disk_end_front_arr, x_disk_end_back_arr, F_blade_arr, material,
                 T1, T_m, p_m, part_number):
        InitialData.__init__(self, n, z_blade, r_arr, h_arr, x_disk_end_front_arr, x_disk_end_back_arr, F_blade_arr,
                             material, T1, T_m, p_m)
        self.part_number = part_number   # число участков разбиения
        self.r1 = self.r_arr_init[0]
        self.r_m = self.r_arr_init[len(self.r_arr_init) - 1]
        self.r_arr = np.linspace(self.r1, self.r_m, part_number + 1)
        self.theta0 = lambda r: 0
        self.theta0_arr = [self.theta0(i) for i in self.r_arr]
        self.E_arr = [self.E(i) for i in self.r_arr]
        self.theta_arr = [self.theta(i) for i in self.r_arr]
        self.h_arr = [self.h(i) for i in self.r_arr]
        self.x_disk_end_front_arr = [self.x_disk_end_front(i) for i in self.r_arr]
        self.x_disk_end_back_arr = [self.x_disk_end_back(i) for i in self.r_arr]
        self.rho_new_arr = [self.rho_new(i) for i in self.r_arr]
        self.sigma_02_arr = [self.sigma_02(i) for i in self.r_arr]

    def rho_new(self, r):
        if r == 0 and self.F_blade(r) == 0:
            return self.material.rho
        elif r == 0 and self.F_blade(r) == 0:
            assert r != 0
        else:
            return self.material.rho * (1 + self.F_blade(r) * 0.5 * self.z_blade / (2 * np.pi * r * self.h(r)))

    def T(self, r):
        return self.T1 + (self.T_m - self.T1) * ((r - self.r1) / (self.r_m - self.r1)) ** 2

    def E(self, r):
        return self.material.E(self.T(r))

    def alpha(self, r):
        return self.material.alpha(self.T(r))

    def theta(self, r):
        return self.alpha(r) * (self.T(r) - 20)

    def sigma_02(self, r):
        return self.material.sigma_02(self.T(r))