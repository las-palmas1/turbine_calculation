import logging
import os
import functions as func
from average_streamline_calculation.turbine import Turbine
from average_streamline_calculation.turbine_stage_geometry import StageGeomAndHeatDrop, InvalidStageSizeValue
from average_streamline_calculation.turbine_stage_gas_dynamics import StageGasDynamics
import enum
import numpy as np
import matplotlib.pyplot as plt


log_filename = os.path.join(os.path.dirname(__file__), 'profiling.log')
logger = func.create_logger(__name__, logging.DEBUG, filename=log_filename, filemode='a',
                            add_console_handler=False)


class ProfilingType(enum.Enum):
    ConstantCirculation = 0
    ConstantAngle = 1


class StageParametersRadialDistribution:
    def __init__(self, profiling_type: ProfilingType, p0_stag, H0, T0_stag, phi, psi, c_p, k, D1_in, D1_av,
                 D1_out, D2_in, D2_av, D2_out, c1_av, alpha1_av, p2_av, c2_a_av, n):
        self.profiling_type = profiling_type
        self.p0_stag = p0_stag
        self.H0 = H0
        self.T0_stag = T0_stag
        self.phi = phi
        self.psi = psi
        self.c_p = c_p
        self.k = k
        self.D1_in = D1_in
        self.D1_av = D1_av
        self.D1_out = D1_out
        self.D2_in = D2_in
        self.D2_av = D2_av
        self.D2_out = D2_out
        self.c1_av = c1_av
        self.alpha1_av = alpha1_av
        self.n = n
        self.p2_av = p2_av
        self.c2_a_av = c2_a_av

    def c1_u(self, r):
        c1_u_av = self.c1_av * np.cos(self.alpha1_av)
        if self.profiling_type == ProfilingType.ConstantCirculation:
            return 0.5 * self.D1_av * c1_u_av / r
        elif self.profiling_type == ProfilingType.ConstantAngle:
            return c1_u_av * (0.5 * self.D1_av / r) ** (np.cos(self.alpha1_av) ** 2)

    def c1_a(self, r):
        c1_a_av = self.c1_av * np.sin(self.alpha1_av)
        if self.profiling_type == ProfilingType.ConstantCirculation:
            return c1_a_av
        elif self.profiling_type == ProfilingType.ConstantAngle:
            return c1_a_av * (0.5 * self.D1_av / r) ** (np.cos(self.alpha1_av) ** 2)

    def p2(self, r):
        return self.p2_av

    def c2_a(self, r):
        return self.c2_a_av

    def c1(self, r):
        return np.sqrt(self.c1_a(r) ** 2 + self.c1_u(r) ** 2)

    def alpha1(self, r):
        if self.c1_a(r) / self.c1(r) > 1:
            raise InvalidStageSizeValue('c1_a must be less than c1')
        return np.arcsin(self.c1_a(r) / self.c1(r))

    def H_s(self, r):
        return self.c1(r) ** 2 / (2 * self.phi)

    def p1(self, r):
        return self.p0_stag * (1 - self.H_s(r) / (self.T0_stag * self.c_p)) ** (self.k / (self.k - 1))

    def T1_ad(self, r):
        return self.T0_stag - self.H_s(r) / self.c_p

    def T1(self, r):
        return self.T0_stag - self.H_s(r) * self.phi ** 2 / self.c_p

    def H_l(self, r):
        return self.c_p * self.T0_stag * (1 - (self.p2(r) / self.p1(r)) ** ((self.k - 1) / self.k))

    def u1(self, r):
        return 2 * np.pi * r * self.n / 60

    def u2(self, r):
        return 2 * np.pi * r * self.n / 60

    def w1(self, r):
        a = self.c1(r) ** 2 + self.u1(r) ** 2 - 2 * self.u1(r) * self.c1(r) * self.alpha1(r)
        if a < 0:
            raise InvalidStageSizeValue('w1 can not be calculated')
        return np.sqrt(self.c1(r) ** 2 + self.u1(r) ** 2 - 2 * self.u1(r) * self.c1(r) * np.cos(self.alpha1(r)))

    def w2(self, r):
        return self.psi * np.sqrt(self.w1(r) ** 2 + 2 * self.H_l(r) + self.u2(r) ** 2 - self.u1(r) ** 2)

    def beta2(self, r):
        if self.c2_a(r) / self.w2(r) > 1:
            raise InvalidStageSizeValue('c2_a must be less than w2')
        return np.arcsin(self.c2_a(r) / self.w2(r))

    def beta1(self, r):
        if self.c1(r) * np.cos(self.alpha1(r)) - self.u1(r) >= 0:
            return np.arctan(self.c1_a(r) / (self.c1(r) * np.cos(self.alpha1(r)) - self.u1(r)))
        else:
            return np.pi + np.arctan(self.c1_a(r) / (self.c1(r) * np.cos(self.alpha1(r)) - self.u1(r)))

    def c2_u(self, r):
        return self.w2(r) * np.cos(self.beta2(r)) - self.u2(r)

    def alpha2(self, r):
        if self.c2_u(r) >= 0:
            return np.arctan(self.c2_a(r) / self.c2_u(r))
        else:
            return np.pi + np.arctan(self.c2_a(r) / self.c2_u(r))

    def c2(self, r):
        return np.sqrt(self.c2_a(r) ** 2 + self.c2_u(r) ** 2)

    def rho(self, r):
        return self.H_l(r) * self.T1_ad(r) / (self.H0 * self.T1(r))

    def plot_parameter_distribution(self, par_name: str, figsize=(9, 7)):
        r_in = 0.5 * self.D1_in
        r_out = 0.5 * self.D1_out
        r_av = 0.5 * self.D1_av
        if par_name.find('1') != -1 or par_name == 'H_s':
            r_in = 0.5 * self.D1_in
            r_out = 0.5 * self.D1_out
            r_av = 0.5 * self.D1_av
        elif par_name.find('2') != -1 or par_name == 'H_l':
            r_in = 0.5 * self.D2_in
            r_out = 0.5 * self.D2_out
            r_av = 0.5 * self.D2_av
        get_atr = object.__getattribute__
        par = get_atr(self, par_name)
        y = np.array(np.linspace(r_in, r_out, 100)) / r_av
        deg = np.pi / 180
        x = [par(i) for i in y * r_av]
        if par_name.find('alpha') != -1 or par_name.find('beta') != -1:
            x = [i / deg for i in x]
        plt.figure(figsize=figsize)
        plt.plot(x, y, linewidth=2, color='blue')
        plt.xlabel(par_name, fontsize=16)
        plt.ylabel(r'$\frac{r}{r_{av}}$', fontsize=22)
        plt.grid()
        plt.show()

    def plot_velocity_triangles(self, r_rel=(0, 0.5, 1), figsize=(8, 8)):
        r_arr = [0.5 * (self.D1_in + i * (self.D1_out - self.D1_in)) for i in r_rel]
        title = [r'$r_{rel} = %s$' % i for i in r_rel]
        for n, i in enumerate(r_arr):
            plt.figure(figsize=figsize)
            x_in = np.array([0, -self.c1_u(i), -self.c1_u(i) + self.u1(i), 0])
            y_in = np.array([self.c1_a(i), 0, 0, self.c1_a(i)])
            x_out = np.array([0, self.c2_u(i), self.c2_u(i) + self.u2(i), 0])
            y_out = np.array([self.c1_a(i), self.c1_a(i) - self.c2_a(i), self.c1_a(i) - self.c2_a(i), self.c1_a(i)])
            plt.plot(x_in, y_in, linewidth=2, color='red', label='inlet')
            plt.plot(x_out, y_out, linewidth=2, color='blue', label='outlet')
            plt.xlim(-self.c1_u(i), self.c2_u(i) + self.u2(i))
            plt.ylim(-max(self.c1_a(i), self.c1_u(i)), max(self.c1_a(i), self.c2_u(i) + self.u2(i)))
            plt.grid()
            plt.title(title[n], fontsize=20)
            plt.legend()
            plt.show()

if __name__ == '__main__':
    rad_dist = StageParametersRadialDistribution(ProfilingType.ConstantCirculation, 482591.6, 64790, 1167.4, 0.97,
                                                 0.97, 1184, 1.32, 0.12579, 0.32100, 0.51621, 0.11446, 0.32848,
                                                 0.542513, 249, 0.8716197, 395993, 183.5, 14687.58)
    rad_dist.plot_parameter_distribution('rho')
    rad_dist.plot_velocity_triangles()









