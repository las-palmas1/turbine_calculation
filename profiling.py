import logging
import os
import functions as func
from average_streamline_calculation.turbine import Turbine
from average_streamline_calculation.turbine_stage_geometry import StageGeomAndHeatDrop, InvalidStageSizeValue
from average_streamline_calculation.turbine_stage_gas_dynamics import StageGasDynamics
import enum
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad


log_filename = os.path.join(os.path.dirname(__file__), 'profiling.log')
logger = func.create_logger(__name__, logging.DEBUG, filename=log_filename, filemode='a',
                            add_console_handler=False)


class ProfilingType(enum.Enum):
    ConstantCirculation = 0
    ConstantAngle = 1


class StageParametersRadialDistribution2:
    def __init__(self, profiling_type: ProfilingType, p0_stag, T0_stag, phi, psi, c_p, k, D1_in, D1_av, D1_out, n,
                 c1_av, alpha1_av, T2_stag_av, L_u_av):
        self.profiling_type = profiling_type
        self.p0_stag = p0_stag
        self.T0_stag = T0_stag
        self.phi = phi
        self.psi = psi
        self.c_p = c_p
        self.k = k
        self.D1_in = D1_in
        self.D1_av = D1_av
        self.D1_out = D1_out
        self.c1_av = c1_av
        self.alpha1_av = alpha1_av
        self.n = n
        self.T2_stag_av = T2_stag_av
        self.L_u_av = L_u_av
        self.R = self.c_p * (self.k - 1) / self.k

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

    def u(self, r):
        return 2 * np.pi * r * self.n / 60

    def L_u(self, r):
        return self.L_u_av

    def c2_u(self, r):
        return (self.L_u(r) - self.c1_u(r) * self.u(r)) / self.u(r)

    def c2_a(self, r):
        def int_func(r):
            return self.c2_u(r) ** 2 / r
        return np.sqrt(self.c1_a(0.5 * self.D1_av) ** 2 + self.c1_u(0.5 * self.D1_av) ** 2 - self.c2_u(r) ** 2 -
                       2 * quad(int_func, 0.5 * self.D1_av, r)[0])

    def c2(self, r):
        return np.sqrt(self.c2_a(r) ** 2 + self.c2_u(r) ** 2)

    def alpha2(self, r):
        if self.c2_u(r) >= 0:
            return np.arctan(self.c2_a(r) / self.c2_u(r))
        else:
            return np.pi + np.arctan(self.c2_a(r) / self.c2_u(r))

    def w2_u(self, r):
        return self.c2_u(r) + self.u(r)

    def w2(self, r):
        return np.sqrt(self.w2_u(r) ** 2 + self.c2_a(r) ** 2)

    def w1(self, r):
        a = self.c1(r) ** 2 + self.u(r) ** 2 - 2 * self.u(r) * self.c1(r) * self.alpha1(r)
        if a < 0:
            raise InvalidStageSizeValue('w1 can not be calculated')
        return np.sqrt(self.c1(r) ** 2 + self.u(r) ** 2 - 2 * self.u(r) * self.c1(r) * np.cos(self.alpha1(r)))

    def beta1(self, r):
        if self.c1(r) * np.cos(self.alpha1(r)) - self.u(r) >= 0:
            return np.arctan(self.c1_a(r) / (self.c1(r) * np.cos(self.alpha1(r)) - self.u(r)))
        else:
            return np.pi + np.arctan(self.c1_a(r) / (self.c1(r) * np.cos(self.alpha1(r)) - self.u(r)))

    def beta2(self, r):
        if self.c2_a(r) / self.w2(r) > 1:
            raise InvalidStageSizeValue('c2_a must be less than w2')
        return np.arcsin(self.c2_a(r) / self.w2(r))

    def T2_stag(self, r):
        return self.T2_stag_av

    def H_l(self, r):
        return 0.5 * ((self.w2(r) / self.psi) ** 2 - self.w1(r) ** 2)

    def p2(self, r):
        return self.p1(r) * (1 - self.H_l(r) / (self.c_p * self.T1(r))) # ** (self.k / (self.k - 1))

    def T2(self, r):
        return self.T2_stag(r) - self.c2(r) ** 2 / (2 * self.c_p)

    def T2_check(self, r):
        return self.T1(r) - (self.w2(r) ** 2 - self.w1(r) ** 2) / (2 * self.c_p)

    def H0(self, r):
        return self.c_p * self.T0_stag * (1 - (self.p0_stag / self.p2(r)) ** ((1 - self.k) / self.k))

    def rho(self, r):
        return self.H_l(r) * self.T1_ad(r) / (self.H0(r) * self.T1(r))

    def M_c1(self, r):
        return self.c1(r) / np.sqrt(self.k * self.R * self.T1(r))

    def M_w2(self, r):
        return self.w2(r) / np.sqrt(self.k * self.R * self.T2(r))

    def plot_parameter_distribution(self, par_name: str, figsize=(9, 7), color='blue'):
        r_in = 0.5 * self.D1_in
        r_out = 0.5 * self.D1_out
        r_av = 0.5 * self.D1_av
        get_atr = object.__getattribute__
        par = get_atr(self, par_name)
        y = np.array(np.linspace(r_in, r_out, 100)) / r_av
        deg = np.pi / 180
        x = [par(i) for i in y * r_av]
        if par_name.find('alpha') != -1 or par_name.find('beta') != -1:
            x = [i / deg for i in x]
        # plt.figure(figsize=figsize)
        plt.plot(x, y, linewidth=2, color=color)
        plt.xlabel(par_name, fontsize=16)
        plt.ylabel(r'$\frac{r}{r_{av}}$', fontsize=22)
        plt.grid()
        # plt.show()

    def plot_velocity_triangles(self, r_rel=(0, 0.5, 1), figsize=(8, 8)):
        r_arr = [0.5 * (self.D1_in + i * (self.D1_out - self.D1_in)) for i in r_rel]
        title = [r'$r_{rel} = %s$' % i for i in r_rel]
        for (n, i) in enumerate(r_arr):
            plt.figure(figsize=figsize)
            x_in = np.array([0, -self.c1_u(i), -self.c1_u(i) + self.u(i), 0])
            y_in = np.array([self.c1_a(i), 0, 0, self.c1_a(i)])
            x_out = np.array([0, self.c2_u(i), self.c2_u(i) + self.u(i), 0])
            y_out = np.array([self.c1_a(i), self.c1_a(i) - self.c2_a(i), self.c1_a(i) - self.c2_a(i), self.c1_a(i)])
            plt.plot(x_in, y_in, linewidth=2, color='red', label='inlet')
            plt.plot(x_out, y_out, linewidth=2, color='blue', label='outlet')
            plt.xlim(-self.c1_u(i), self.c2_u(i) + self.u(i))
            plt.ylim(-max(self.c1_a(i), self.c1_u(i)), max(self.c1_a(i), self.c2_u(i) + self.u(i)))
            plt.grid()
            plt.title(title[n], fontsize=20)
            plt.legend()
            plt.show()


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
        self.R = self.c_p * (self.k - 1) / self.k

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
        return self.psi * np.sqrt(self.w1(r) ** 2 + 2 * self.H_l(r) + self.u2(r) ** 2 -
                                  self.u1(r) ** 2)

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

    def T2(self, r):
        return self.T1(r) - (self.w2(r)**2 - self.w1(r)**2 - self.u2(r)**2 + self.u1(r)**2) / (2 * self.c_p)

    def M_c1(self, r):
        return self.c1(r) / np.sqrt(self.k * self.R * self.T1(r))

    def M_w2(self, r):
        return self.w2(r) / np.sqrt(self.k * self.R * self.T2(r))

    def T2_stag(self, r):
        return self.T2(r) + self.c2(r) ** 2 / (2 * self.c_p)

    def L_u(self, r):
        return self.u1(r) * self.c1_u(r) + self.u2(r) * self.c2_u(r)

    def plot_parameter_distribution(self, par_name: str, figsize=(9, 7), color='blue'):
        r_in = 0.5 * self.D1_in
        r_out = 0.5 * self.D1_out
        r_av = 0.5 * self.D1_av
        get_atr = object.__getattribute__
        par = get_atr(self, par_name)
        y = np.array(np.linspace(r_in, r_out, 100)) / r_av
        deg = np.pi / 180
        x = [par(i) for i in y * r_av]
        if par_name.find('alpha') != -1 or par_name.find('beta') != -1:
            x = [i / deg for i in x]
        # plt.figure(figsize=figsize)
        plt.plot(x, y, linewidth=2, color=color)
        plt.xlabel(par_name, fontsize=16)
        plt.ylabel(r'$\frac{r}{r_{av}}$', fontsize=22)
        plt.grid()
        # plt.show()

    def plot_velocity_triangles(self, r_rel=(0, 0.5, 1), figsize=(8, 8)):
        r_arr = [0.5 * (self.D1_in + i * (self.D1_out - self.D1_in)) for i in r_rel]
        title = [r'$r_{rel} = %s$' % i for i in r_rel]
        for (n, i) in enumerate(r_arr):
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


class BladeProfile:
    def __init__(self, angle1, angle2, b_a, gamma1_s, gamma1_k, pnt_count: int, r1, s):
        """

        :param angle1: угол лопатки на входе
        :param angle2: угол лопакти на выходе
        :param b_a:
        :param gamma1_s:
        :param gamma1_k:
        :param pnt_count:
        :param r1: рвдиус скругления входной кромки
        :param s: толщина выходной кромки
        """
        self._angle1 = angle1
        self._angle2 = angle2
        self._b_a = b_a
        self._gamma1_s = gamma1_s
        self._gamma1_k = gamma1_k
        self._pnt_count = pnt_count
        self._r1 = r1
        self._s = s
        self._compute_coordinates()

    @property
    def s(self) -> int:
        assert self._s is not None, 's must not be None'
        return self._s

    @s.setter
    def s(self, value: int):
        self._s = value
        self._compute_coordinates()

    @property
    def r1(self) -> int:
        assert self._r1 is not None, 'r1 must not be None'
        return self._r1

    @r1.setter
    def r1(self, value: int):
        self._r1 = value
        self._compute_coordinates()

    @property
    def pnt_count(self) -> int:
        assert self._pnt_count is not None, 'pnt_count must not be None'
        return self._pnt_count

    @pnt_count.setter
    def pnt_count(self, value: int):
        self._pnt_count = value
        self._compute_coordinates()

    @property
    def gamma1_k(self):
        assert self._gamma1_k is not None, 'gamma1_k must not be None'
        return self._gamma1_k

    @gamma1_k.setter
    def gamma1_k(self, value):
        self._gamma1_k = value
        self._compute_coordinates()

    @property
    def gamma1_s(self):
        assert self._gamma1_s is not None, 'gamma1_s must not be None'
        return self._gamma1_s

    @gamma1_s.setter
    def gamma1_s(self, value):
        self._gamma1_s = value
        self._compute_coordinates()

    @property
    def b_a(self):
        assert self._b_a is not None, 'b_a must not be None'
        return self._b_a

    @b_a.setter
    def b_a(self, value):
        self._b_a = value
        self._compute_coordinates()

    @property
    def angle2(self):
        assert self._angle2 is not None, 'angle2 must not be None'
        return self._angle2

    @angle2.setter
    def angle2(self, value):
        self._angle2 = value
        self._compute_coordinates()

    @property
    def angle1(self):
        assert self._angle1 is not None, 'angle1 must not be None'
        return self._angle1

    @angle1.setter
    def angle1(self, value):
        self._angle1 = value
        self._compute_coordinates()

    @classmethod
    def compute_parabola_coordinates_by_dir(cls, x1, x2, dir1, dir2, pnt_count: int, c=0):
        a = (dir2 - dir1) / (2 * (x2 - x1))
        b = (x2 * dir1 - x1 * dir2) / (x2 - x1)
        x = np.array(np.linspace(x1, x2, pnt_count))
        y = a * x ** 2 + b * x + c
        return x, y

    @classmethod
    def compute_parabola_coordinates_by_points(cls, x1, y1, x2, y2, pnt_count, **kwargs):
        a, b, c = None, None, None
        if 'dir1' in kwargs:
            A = [[x1 ** 2, x1, 1], [x2 ** 2, x2, 1], [2 * x1, 1, 0]]
            B = [y1, y2, kwargs['dir1']]
            a, b, c = la.solve(A, B)
        elif 'dir2' in kwargs:
            A = [[x1 ** 2, x1, 1], [x2 ** 2, x2, 1], [2 * x2, 1, 0]]
            B = [y1, y2, kwargs['dir2']]
            a, b, c = la.solve(A, B)
        x = np.array(np.linspace(x1, x2, pnt_count))
        y = a * x ** 2 + b * x + c
        return x, y

    def _compute_coordinates(self):
        try:
            self.angle1_s = self.angle1 - self.gamma1_s
            self.angle1_k = self.angle1 + self.gamma1_k
            self.dir1_s = np.tan(self.angle1_s + np.pi / 2)
            self.dir1_k = np.tan(self.angle1_k + np.pi / 2)
            self.dir1_av = np.tan(self.angle1 + np.pi / 2)
            self.dir2_av = np.tan(np.pi / 2 - self.angle2)
            self.x_av, self.y_av = self.compute_parabola_coordinates_by_dir(0, self.b_a, self.dir1_av,
                                                                            self.dir2_av, self.pnt_count)
            self.x_s, self.y_s = self.compute_parabola_coordinates_by_points(0, self.y_av[0] - 0.5 * self.r1, self.b_a,
                                                                             self.y_av[self.pnt_count - 1] -
                                                                             0.5 * self.s, self.pnt_count,
                                                                             dir1=self.dir1_s)
            self.x_k, self.y_k = self.compute_parabola_coordinates_by_points(0, self.y_av[0] + 0.5 * self.r1, self.b_a,
                                                                             self.y_av[self.pnt_count - 1] +
                                                                             0.5 * self.s, self.pnt_count,
                                                                             dir1=self.dir1_k)
            self.angle2_s = -np.arctan((self.y_s[self.pnt_count - 1] - self.y_s[self.pnt_count - 2]) /
                                       (self.x_s[self.pnt_count - 1] - self.x_s[self.pnt_count - 2])) + np.pi / 2
            self.angle2_k = -np.arctan((self.y_k[self.pnt_count - 1] - self.y_k[self.pnt_count - 2]) /
                                       (self.x_k[self.pnt_count - 1] - self.x_k[self.pnt_count - 2])) + np.pi / 2
            self.gamma2_s = self.angle2 - self.angle2_s
            self.gamma2_k = self.angle2_k - self.angle2
        except AssertionError:
            pass


class StageProfiling(StageParametersRadialDistribution):
    def __init__(self, profiling_type: ProfilingType, p0_stag, H0, T0_stag, phi, psi, c_p, k, D1_in, D1_av,
                 D1_out, D2_in, D2_av, D2_out, c1_av, alpha1_av, p2_av, c2_a_av, n, alpha0, pnt_count, b_a_sa,
                 b_a_rk, r_rel=(0, 0.5, 1), gamma1_k_sa_rel=0.25, gamma1_k_rk_rel=0.25, s2_sa=0.001, s2_rk=0.001):
        StageParametersRadialDistribution.__init__(self, profiling_type, p0_stag, H0, T0_stag, phi, psi, c_p, k, D1_in,
                                                   D1_av, D1_out, D2_in, D2_av, D2_out, c1_av, alpha1_av, p2_av,
                                                   c2_a_av, n)
        self.alpha0 = alpha0
        self.r_rel = r_rel
        self.pnt_count = pnt_count
        self.b_a_sa = b_a_sa
        self.b_a_rk = b_a_rk
        self.gamma1_k_sa_rel = gamma1_k_sa_rel
        self.gamma1_k_rk_rel = gamma1_k_rk_rel
        self.s2_sa = s2_sa
        self.s2_rk = s2_rk
        self.r1_sa = 0.04 * self.b_a_sa
        self.r1_rk = 0.04 * self.b_a_rk
        self.gamma_sa = None
        self.gamma_rk = None
        self.t_rel_sa = None
        self.t_rel_rk = None
        self.b_sa = None
        self.b_rk = None
        self.z_sa = None
        self.z_rk = None
        self.t_rk_av = None
        self.t_sa_av = None

    @classmethod
    def delta(cls, M, angle1):
        """
        :param M: число Маха
        :return:
        Возвращает значение угла оставания
        """
        deg = np.pi / 180
        if M >= 0.95:
            return 0
        elif 15 * deg <= angle1 <= 20 * deg:
            M = [0.95, 0.8]
            delta = [1 / 3 * deg, 5 / 6 * deg]
            delta_int = interp1d(M, delta)
            return delta_int(M)
        elif 20 * deg < angle1 <= 30 * deg:
            M = [0.95, 0.6]
            delta = [1 / 3 * deg, 1.5 * deg]
            delta_int = interp1d(M, delta)
            return delta_int(M)
        elif 30 * deg < angle1 <= 40 * deg:
            M = [0.95, 0.6]
            delta = [1 / 3 * deg, 4 * deg]
            delta_int = interp1d(M, delta)
            return delta_int(M)
        else:
            return 1 * deg

    @classmethod
    def gamma1(cls, angle1):
        deg = np.pi / 180
        angle1_arr = np.array([10, 30, 150]) * deg
        gamma1_arr = np.array([55, 45, 10]) * deg
        gamma1_int = interp1d(angle1_arr, gamma1_arr)
        return gamma1_int(angle1)

    def alpha0_l(self, r):
        return self.alpha0(r)

    def alpha1_l(self, r):
        return self.alpha1(r) - self.delta(self.M_c1(r), self.alpha1(r))

    def beta1_l(self, r):
        return self.beta1(r)

    def beta2_l(self, r):
        return self.beta2(r) - self.delta(self.M_w2(r), self.beta2(r))

    def compute_step(self):
        deg = np.pi / 180
        self.gamma_rk = (70 - 0.127 * (self.beta1(0.5 * self.D2_av) - self.beta2(0.5 * self.D2_av)) / deg -
                         0.0041 * (self.beta1(0.5 * self.D2_av) - self.beta2(0.5 * self.D2_av)) ** 2 / deg ** 2) * deg
        self.gamma_sa = (70 - 0.127 * (self.alpha0(0.5 * self.D1_av) - self.alpha1(0.5 * self.D1_av)) / deg -
                         0.0041 * (self.alpha0(0.5 * self.D1_av) - self.alpha1(0.5 * self.D1_av)) ** 2 / deg ** 2) * deg
        self.b_rk = self.b_a_rk / np.sin(self.gamma_rk)
        self.b_sa = self.b_a_sa / np.sin(self.gamma_sa)
        self.t_rel_rk = 0.6 * (180 / (180 - (self.beta1(0.5 * self.D2_av) + self.beta2(0.5 * self.D2_av)) / deg) *
                               np.sin(self.beta1(0.5 * self.D2_av)) / np.sin(self.beta2(0.5 * self.D2_av))) ** \
                              (1 / 3) * (1 - 0.13)
        self.t_rel_sa = 0.6 * (180 / (180 - (self.alpha0(0.5 * self.D1_av) + self.alpha1(0.5 * self.D1_av)) / deg) *
                               np.sin(self.alpha0(0.5 * self.D1_av)) / np.sin(self.alpha1(0.5 * self.D1_av))) ** \
                              (1 / 3) * (1 - 0.15)
        self.z_rk = int(np.pi * self.D2_av / (self.t_rel_rk * self.b_rk))
        self.z_sa = int(np.pi * self.D1_av / (self.t_rel_sa * self.b_sa))
        self.t_rk_av = np.pi * self.D2_av / self.z_rk
        self.t_sa_av = np.pi * self.D1_av / self.z_sa

    def gamma1_k_sa(self, r):
        return self.gamma1(self.alpha0(r)) * self.gamma1_k_sa_rel

    def gamma1_s_sa(self, r):
        return self.gamma1(self.alpha0(r)) * (1 - self.gamma1_k_sa_rel)

    def gamma1_k_rk(self, r):
        return self.gamma1(self.beta1(r)) * self.gamma1_k_rk_rel

    def gamma1_s_rk(self, r):
        return self.gamma1(self.beta1(r)) * (1 - self.gamma1_k_rk_rel)


if __name__ == '__main__':
    rad_dist1 = StageParametersRadialDistribution(ProfilingType.ConstantCirculation, 482591.6, 64790, 1167.4, 0.97,
                                                  0.97, 1184, 1.32, 0.12579, 0.32100, 0.51621, 0.11446, 0.32848,
                                                  0.542513, 249, 0.8716197, 395993, 183.5, 14687.58)
    rad_dist2 = StageParametersRadialDistribution(ProfilingType.ConstantAngle, 482591.6, 64790, 1167.4, 0.97,
                                                  0.97, 1184, 1.32, 0.12579, 0.32100, 0.51621, 0.11446, 0.32848,
                                                  0.542513, 249, 0.8716197, 395993, 183.5, 14687.58)

    plt.figure()
    rad_dist1.plot_parameter_distribution('rho', color='red')
    rad_dist2.plot_parameter_distribution('rho')
    plt.xlim()
    plt.show()

    rad_dist3 = StageParametersRadialDistribution2(ProfilingType.ConstantCirculation, 482591.6, 1167.4, 0.97, 0.97,
                                                   1184, 1.32, 0.12579, 0.32100, 0.51621, 14687.58, 249, 0.8716197,
                                                   1132.57, 42704.509)
    rad_dist4 = StageParametersRadialDistribution2(ProfilingType.ConstantAngle, 482591.6, 1167.4, 0.97, 0.97,
                                                   1184, 1.32, 0.12579, 0.32100, 0.51621, 14687.58, 249, 0.8716197,
                                                   1132.57, 42704.509)
    plt.figure()
    rad_dist3.plot_parameter_distribution('c2_u', color='red')
    rad_dist4.plot_parameter_distribution('c2_u')
    plt.xlim()
    plt.show()
    print(rad_dist3.c2_u(rad_dist3.D1_av / 2))
    # deg = np.pi / 180
    # bl_prof = BladeProfile(45 * deg, 20 * deg, 0.06, 20 * deg, 15 * deg, 100, 0.003, 0.001)
    # plt.figure(figsize=(8, 8))
    # plt.plot(bl_prof.y_av, -bl_prof.x_av, color='black')
    # plt.plot(bl_prof.y_s, -bl_prof.x_s, color='red',)
    # plt.plot(bl_prof.y_k, -bl_prof.x_k, color='green')
    # plt.grid()
    # plt.show()
    # print(bl_prof.gamma2_s / deg)
    # print(bl_prof.gamma2_k / deg)








