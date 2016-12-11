import logging
import os
import functions as func
from gas_dynamics import GasDynamicFunctions
from average_streamline_calculation.turbine import Turbine, TurbineType
from average_streamline_calculation.turbine_stage_geometry import StageGeomAndHeatDrop, InvalidStageSizeValue
from average_streamline_calculation.turbine_stage_gas_dynamics import StageGasDynamics
import enum
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad
from mpl_toolkits.mplot3d import Axes3D
import pickle as pk


log_filename = os.path.join(os.path.dirname(__file__), 'profiling.log')
logger = func.create_logger(__name__, logging.DEBUG, filename=log_filename, filemode='a',
                            add_console_handler=False)


class ProfilingType(enum.Enum):
    ConstantCirculation = 0
    ConstantAngle = 1


class StageParametersRadialDistribution:
    def __init__(self, profiling_type: ProfilingType, p0_stag, T0_stag, c0, c_p, k, D1_in, D1_av, D1_out, n,
                 c1_av, alpha1_av, T2_stag_av, L_u_av, c2_a_av, c2_u_av, phi=1.0, psi=1.0):
        self.profiling_type = profiling_type
        self.p0_stag = p0_stag
        self.T0_stag = T0_stag
        self.c0 = c0
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
        self.c2_a_av = c2_a_av
        self.c2_u_av = c2_u_av
        self.R = self.c_p * (self.k - 1) / self.k

    def T0(self, r):
        return self.T0_stag - self.c0(r) ** 2 / (2 * self.c_p)

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
        return self.c1(r) ** 2 / (2 * self.phi ** 2)

    def p1(self, r):
        return self.p0_stag(r) * (1 - self.H_s(r) / (self.T0_stag * self.c_p)) ** (self.k / (self.k - 1))

    def T1_ad(self, r):
        return self.T0_stag - self.H_s(r) / self.c_p

    def T1(self, r):
        return self.T0_stag - self.H_s(r) * self.phi ** 2 / self.c_p

    def p1_ideal(self, r):
        def int_func(r):
            return self.c1_u(r) ** 2 / (self.R * self.T1(r) * r)
        return self.p1(self.D1_av / 2) * np.exp(quad(int_func, 0.5 * self.D1_av, r)[0])

    def u(self, r):
        return 2 * np.pi * r * self.n / 60

    def L_u(self, r):
        return self.L_u_av

    def c2_u(self, r):
        return (self.L_u(r) - self.c1_u(r) * self.u(r)) / self.u(r)

    def c2_a(self, r):
        def int_func(r):
            return self.c2_u(r) ** 2 / r

        return np.sqrt(self.c2_a_av ** 2 + self.c2_u_av ** 2 - self.c2_u(r) ** 2 -
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
        return self.p1(r) * (1 - self.H_l(r) / (self.c_p * self.T1(r))) ** (self.k / (self.k - 1))

    def p2_ideal(self, r):
        def int_func(r):
            return self.c2_u(r) ** 2 / (self.R * self.T2(r) * r)
        return self.p2(self.D1_av / 2) * np.exp(quad(int_func, 0.5 * self.D1_av, r)[0])

    def p2_stag(self, r):
        return self.p2(r) / GasDynamicFunctions.pi_lam(self.c2(r) /
                                                       GasDynamicFunctions.a_cr(self.T2_stag(r), self.k, self.R), self.k)

    def T2(self, r):
        return self.T1(r) - (self.w2(r) ** 2 - self.w1(r) ** 2) / (2 * self.c_p)

    def H0(self, r):
        return self.c_p * self.T0_stag * (1 - (self.p0_stag(r) / self.p2(r)) ** ((1 - self.k) / self.k))

    def rho(self, r):
        return self.H_l(r) * self.T1_ad(r) / (self.H0(r) * self.T1(r))

    def M_c1(self, r):
        return self.c1(r) / np.sqrt(self.k * self.R * self.T1(r))

    def M_w2(self, r):
        return self.w2(r) / np.sqrt(self.k * self.R * self.T2(r))

    def plot_parameter_distribution(self, par_names, color, figsize=(9, 7)):
        r_in = 0.5 * self.D1_in
        r_out = 0.5 * self.D1_out
        r_av = 0.5 * self.D1_av
        get_atr = object.__getattribute__
        y1 = np.array(np.linspace(r_in, r_out, 100)) / r_av
        y = np.array(np.linspace(r_in, r_out, 100))
        deg = np.pi / 180
        plt.figure(figsize=figsize)
        for n, item in enumerate(par_names):
            par = get_atr(self, item)
            x = [par(i) for i in y]
            if item.find('alpha') != -1 or item.find('beta') != -1 or item.find('delta') != -1 or \
                           item.find('gamma') != -1:
                x = [i / deg for i in x]
            plt.plot(x, y1, linewidth=2, color=color[n], label=item)
        plt.legend(fontsize=16)
        plt.ylabel(r'$\frac{r}{r_{av}}$', fontsize=22)
        plt.grid()
        plt.show()

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


class BladeSection:
    def __init__(self, angle1=None, angle2=None, b_a=None, gamma1_s=None, gamma1_k=None, pnt_count=20, r1=None,
                 s2=0.001):
        """

        :param angle1: угол лопатки на входе
        :param angle2: угол лопакти на выходе
        :param b_a:
        :param gamma1_s:
        :param gamma1_k:
        :param pnt_count:
        :param r1: рвдиус скругления входной кромки
        :param s2: толщина выходной кромки
        """
        self._angle1 = angle1
        self._angle2 = angle2
        self._b_a = b_a
        self._gamma1_s = gamma1_s
        self._gamma1_k = gamma1_k
        self._pnt_count = pnt_count
        self._r1 = r1
        self._s2 = s2
        if (self._angle1 is not None) and (self._angle2 is not None) and (self._b_a is not None) and \
                (self._gamma1_k is not None) and (self._gamma1_s is not None) and (self._pnt_count is not None) and \
                (self._r1 is not None) and (self._s2 is not None):
            self._compute_coordinates()
        self.r = None

    @property
    def s2(self) -> int:
        assert self._s2 is not None, 's2 must not be None'
        return self._s2

    @s2.setter
    def s2(self, value: int):
        self._s2 = value
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

    def _compute_section_center(self):
        x_av_arr, y_av_arr = self.compute_parabola_coordinates_by_dir(0, self.b_a, self.dir1_av, self.dir2_av, 50)
        x_s_arr, y_s_arr = self.compute_parabola_coordinates_by_points(0, y_av_arr[0] - 0.5 * self.r1, self.b_a,
                                                                       y_av_arr[49] - 0.5 * self.s2, 50,
                                                                       dir1=self.dir1_s)
        x_k_arr, y_k_arr = self.compute_parabola_coordinates_by_points(0, y_av_arr[0] + 0.5 * self.r1, self.b_a,
                                                                       y_av_arr[49] + 0.5 * self.s2, 50,
                                                                       dir1=self.dir1_k)

        y_s_int = interp1d(x_s_arr, y_s_arr)
        y_k_int = interp1d(x_k_arr, y_k_arr)

        def y_s(x):
            return y_s_int(x)

        def y_k(x):
            return y_k_int(x)

        self.square = quad(y_k, min(self.x_k), max(self.x_k))[0] - quad(y_s, min(self.x_s), max(self.x_s))[0]

        def func_to_int_sx(x):
            return y_k(x) ** 2 - y_s(x) ** 2

        self.s_x = quad(func_to_int_sx, min(self.x_s), max(self.x_s))[0] * 0.5
        self.y_center = self.s_x / self.square

        x_s_min = x_s_arr[list(y_s_arr).index(min(y_s_arr))]
        x_k_min = x_k_arr[list(y_k_arr).index(min(y_k_arr))]
        x_s_arr1 = np.linspace(min(x_s_arr), x_s_min, 30)
        x_k_arr1 = np.linspace(min(x_k_arr), x_k_min, 30)
        x_s_arr2 = np.linspace(x_s_min, max(x_s_arr), 30)
        x_k_arr2 = np.linspace(x_k_min, max(x_k_arr), 30)
        y_s_arr1 = y_s(x_s_arr1)
        y_k_arr1 = y_k(x_k_arr1)
        y_s_arr2 = y_s(x_s_arr2)
        y_k_arr2 = y_k(x_k_arr2)
        x_s1_int = interp1d(y_s_arr1, x_s_arr1)
        x_s2_int = interp1d(y_s_arr2, x_s_arr2)
        x_k1_int = interp1d(y_k_arr1, x_k_arr1)
        x_k2_int = interp1d(y_k_arr2, x_k_arr2)

        def x_s1(y):
            return x_s1_int(y) ** 2 / 2

        def x_s2(y):
            return x_s2_int(y) ** 2 / 2

        def x_k1(y):
            return x_k1_int(y) ** 2 / 2

        def x_k2(y):
            return x_k2_int(y) ** 2 / 2

        s_y1 = quad(x_k2, y_k_arr[len(y_k_arr) - 1], min(y_k_arr))[0]
        s_y2 = quad(x_k1, min(y_k_arr), y_k_arr[0])[0]
        s_y3 = x_s_arr[0] ** 2 / 2 * (y_s_arr[0] - y_k_arr[0])
        s_y4 = quad(x_s1, y_s_arr[0], min(y_s_arr))[0]
        s_y5 = quad(x_s2, min(y_s_arr), y_s_arr[len(y_s_arr) - 1])[0]
        s_y6 = x_s_arr[len(x_s_arr) - 1] ** 2 / 2 * (y_k_arr[len(y_k_arr) - 1] - y_s_arr[len(y_k_arr) - 1])

        self.s_y = s_y1 + s_y2 + s_y3 + s_y4 + s_y5 + s_y6
        self.x_center = self.s_y / self.square

    def _compute_coordinates(self):
        try:
            self.angle1_s = self.angle1 - self.gamma1_s
            self.angle1_k = self.angle1 + self.gamma1_k
            self.dir1_s = np.tan(self.angle1_s + np.pi / 2)
            self.dir1_k = np.tan(self.angle1_k + np.pi / 2)
            self.dir1_av = np.tan(self.angle1 + np.pi / 2)
            self.dir2_av = np.tan(np.pi / 2 - self.angle2)
            self.delta_y1 = self.r1 * (1 / np.tan(self.angle1_k / 2) + 1 / np.tan(self.angle1_s / 2))
            self.x_av, self.y_av = self.compute_parabola_coordinates_by_dir(0, self.b_a, self.dir1_av,
                                                                            self.dir2_av, self.pnt_count)
            self.x_s, self.y_s = self.compute_parabola_coordinates_by_points(0, self.y_av[0] -
                                                                             0.5 * self.delta_y1,
                                                                             self.b_a, self.y_av[self.pnt_count - 1] -
                                                                             0.5 * self.s2,
                                                                             self.pnt_count, dir1=self.dir1_s)
            self.x_k, self.y_k = self.compute_parabola_coordinates_by_points(0, self.y_av[0] +
                                                                             0.5 * self.delta_y1,
                                                                             self.b_a, self.y_av[self.pnt_count - 1] +
                                                                             0.5 * self.s2,
                                                                             self.pnt_count, dir1=self.dir1_k)
            self._compute_section_center()
            self.x_s -= self.x_center
            self.y_s -= self.y_center
            self.x_k -= self.x_center
            self.y_k -= self.y_center
            self.x_center = 0
            self.y_center = 0
            self.angle2_s = -np.arctan((self.y_s[self.pnt_count - 1] - self.y_s[self.pnt_count - 2]) /
                                       (self.x_s[self.pnt_count - 1] - self.x_s[self.pnt_count - 2])) + np.pi / 2
            self.angle2_k = -np.arctan((self.y_k[self.pnt_count - 1] - self.y_k[self.pnt_count - 2]) /
                                       (self.x_k[self.pnt_count - 1] - self.x_k[self.pnt_count - 2])) + np.pi / 2
            self.gamma2_s = self.angle2 - self.angle2_s
            self.gamma2_k = self.angle2_k - self.angle2
            self.tan_alpha = self.b_a / (self.y_k[len(self.y_av) - 1] - self.y_k[0])
            # alpha - угол установки
            if self.tan_alpha >= 0:
                self.alpha = np.arctan(self.tan_alpha)
            else:
                self.alpha = np.pi + np.arctan(self.tan_alpha)
        except AssertionError:
            pass


class StageProfiling(StageParametersRadialDistribution):
    deg = np.pi / 180

    def __init__(self, profiling_type: ProfilingType, p0_stag, T0_stag, c0, c_p, k, D1_in, D1_av, D1_out, n,
                 c1_av, alpha1_av, T2_stag_av, L_u_av, c2_a_av, c2_u_av, alpha0, pnt_count, b_a_sa, b_a_rk, delta_a_sa,
                 psi=1.0, phi=1.0, r_rel=(0, 0.5, 1), gamma1_k_sa_rel=0.25, gamma1_k_rk_rel=0.25, s2_sa=0.001,
                 s2_rk=0.001, mindif_alpha0_l_gamma1_s_sa=20 * deg, mindif_beta1_l_gamma1_s_rk=17*deg,
                 maxdif_alpha0_l_gamma1_s_sa=40*deg, maxdif_beta1_l_gamma1_s_rk=37*deg, min_t_rel_sa_in=0.75,
                 min_t_rel_rk_in=0.75, r1_rk_rel=0.04, r1_sa_rel=0.04):
        """

        :param profiling_type:
        :param p0_stag:
        :param T0_stag:
        :param phi:
        :param psi:
        :param c_p:
        :param k:
        :param D1_in:
        :param D1_av:
        :param D1_out:
        :param n:
        :param c1_av:
        :param alpha1_av:
        :param T2_stag_av:
        :param L_u_av:
        :param c2_a_av:
        :param c2_u_av:
        :param alpha0: угол входа в СА в зависимости от радиуса
        :param pnt_count:  количесиство точек для нахождения координат сторон профиля в каждом сечении
        :param b_a_sa:
        :param b_a_rk:
        :param r_rel: относительные радиусы сечений для нахождения координат
        :param gamma1_k_sa_rel: отношение углов gamma1_k и gamma1 для СА
        :param gamma1_k_rk_rel: отношение углов gamma1_k и gamma1 для РК
        :param s2_sa: толщина выходной кромки СА
        :param s2_rk: толщина выходной кромки РК
        :param r1_rk_rel: отношение радиуса закругления входной кромки РК к осевой ширине
        :param r1_sa_rel: отношение радиуса закругления входной кромки СА к осевой ширине
        :param mindif_alpha0_l_gamma1_s_sa: минимальная разность углов alpha0_l и gamma1_s_sa
        :param mindif_beta1_l_gamma1_s_rk: минимальная разность углов beta1_l и gamma1_s_rk
        """
        StageParametersRadialDistribution.__init__(self, profiling_type, p0_stag, T0_stag, c0, c_p, k, D1_in,
                                                   D1_av, D1_out, n, c1_av, alpha1_av, T2_stag_av, L_u_av, c2_a_av,
                                                   c2_u_av, phi=phi, psi=psi)
        self.alpha0 = alpha0
        self.r_rel = r_rel
        self.pnt_count = pnt_count
        self.b_a_sa = b_a_sa
        self.b_a_rk = b_a_rk
        self.gamma1_k_sa_rel = gamma1_k_sa_rel
        self.gamma1_k_rk_rel = gamma1_k_rk_rel
        self.s2_sa = s2_sa
        self.s2_rk = s2_rk
        self.r1_sa_rel = r1_sa_rel
        self.r1_rk_rel = r1_rk_rel
        self.mindif_alpha0_l_gamma1_s_sa = mindif_alpha0_l_gamma1_s_sa
        self.mindif_beta1_l_gamma1_s_rk = mindif_beta1_l_gamma1_s_rk
        self.maxdif_alpha0_l_gamma1_s_sa = maxdif_alpha0_l_gamma1_s_sa
        self.maxdif_beta1_l_gamma1_s_rk = maxdif_beta1_l_gamma1_s_rk
        self.min_t_rel_sa_in = min_t_rel_sa_in
        self.min_t_rel_rk_in = min_t_rel_rk_in
        self.delta_a_sa = delta_a_sa
        self.gamma_sa = None
        self.gamma_rk = None
        self.z_sa = None
        self.z_rk = None
        self.t_rk_av = None
        self.t_sa_av = None
        self._sa_sections = [BladeSection(r1=self.r1_sa_rel * self.b_a_sa, s2=self.s2_sa, pnt_count=self.pnt_count,
                                          b_a=self.b_a_sa)
                             for _ in self.r_rel]
        self._rk_sections = [BladeSection(r1=self.r1_rk_rel * self.b_a_rk, s2=self.s2_rk, pnt_count=self.pnt_count,
                                          b_a=self.b_a_rk)
                             for _ in self.r_rel]

    def __getitem__(self, item):
        if 0 <= item < len(self.r_rel):
            return {'sa': self._sa_sections[item], 'rk': self._rk_sections[item]}
        else:
            raise IndexError('invalid index')

    def __len__(self):
        return len(self.r_rel)

    def __iter__(self):
        self._num = 0
        return self

    def __next__(self):
        if self._num < len(self.r_rel):
            current = {'sa': self._sa_sections[self._num], 'rk': self._rk_sections[self._num]}
            self._num += 1
            return current
        else:
            raise StopIteration()

    @classmethod
    def delta_f(cls, M, angle2):
        """
        :param M: число Маха
        :return:
        Возвращает значение угла оставания
        """
        deg = np.pi / 180
        if M >= 0.95:
            return 0
        elif 15 * deg <= angle2 <= 20 * deg:
            M_arr = np.array([0.95, 0.8, 0.3, 0.1])
            delta = np.array([1 / 3 * deg, 5 / 6 * deg, 8 / 6 * deg, 10 / 6 * deg])
            delta_int = interp1d(M_arr, delta)
            return delta_int(M)
        elif 20 * deg < angle2 <= 30 * deg:
            M_arr = np.array([0.95, 0.6, 0.3, 0.1])
            delta = np.array([1 / 3 * deg, 1.5 * deg, 2 * deg, 2.5 * deg])
            delta_int = interp1d(M_arr, delta)
            return delta_int(M)
        elif 30 * deg < angle2 <= 40 * deg:
            M_arr = np.array([0.95, 0.6, 0.2])
            delta = np.array([1 / 3 * deg, 4 * deg, 6 * deg])
            delta_int = interp1d(M_arr, delta)
            return delta_int(M)
        else:
            return 1 * deg

    @classmethod
    def gamma1(cls, angle1):
        deg = np.pi / 180
        angle1_arr = np.array([5, 170]) * deg
        gamma1_arr = np.array([22, 5]) * deg
        gamma1_int = interp1d(angle1_arr, gamma1_arr)
        return gamma1_int(angle1)

    def delta_sa(self, r):
        return self.delta_f(self.M_c1(0.5 * self.D1_av), self.alpha1(0.5 * self.D1_av))

    def delta_rk(self, r):
        return self.delta_f(self.M_w2(0.5 * self.D1_av), self.beta2(0.5 * self.D1_av))

    def alpha0_l(self, r):
        return self.alpha0(r)

    def alpha1_l(self, r):
        return self.alpha1(r) - self.delta_sa(r)

    def beta1_l(self, r):
        return self.beta1(r)

    def beta2_l(self, r):
        return self.beta2(r) - self.delta_rk(r)

    def gamma(self, angle1, angle2):
        deg = np.pi / 180
        return (70 - 0.127 * (angle1 - angle2) / deg - 0.0041 * (angle1 - angle2) ** 2 / deg ** 2) * deg

    def b_sa(self, r):
        return self.b_a_sa / np.sin(self.gamma(self.alpha0(r), self.alpha1(r)))

    def b_rk(self, r):
        return self.b_a_rk / np.sin(self.gamma(self.beta1(r), self.beta2(r)))

    def t_rel_rk_func(self, r):
        return 0.6 * (180 / (180 - (self.beta1(r) + self.beta2(r)) / self.deg) *
                      np.sin(self.beta1(r)) / np.sin(self.beta2(r))) ** (1 / 3) * (1 - 0.11)

    def t_rel_sa_func(self, r):
        return 0.45 * (180 / (180 - (self.alpha0(r) + self.alpha1(r)) / self.deg) *
                       np.sin(self.alpha0(r)) / np.sin(self.alpha1(r))) ** (1 / 3) * (1 - 0.13)

    def _compute_step(self):
        t_rel_rk_av = self.t_rel_rk_func(0.5 * self.D1_av)
        t_rel_sa_av = self.t_rel_sa_func(0.5 * self.D1_av)
        if self.rho(0.5 * self.D1_av) > 0.1:
            if t_rel_rk_av < 0.8:
                t_rel_rk_av = 0.8
            elif t_rel_rk_av > 1.0:
                t_rel_rk_av = 1.0
        elif self.rho(0.5 * self.D1_av) < 0.1:
            if t_rel_rk_av < 0.65:
                t_rel_rk_av = 0.65
            elif t_rel_rk_av > 0.85:
                t_rel_rk_av = 0.85
        if t_rel_sa_av < 0.8:
            t_rel_sa_av = 0.8
        elif t_rel_sa_av > 1.0:
            t_rel_sa_av = 1.0
        self.z_rk = int(np.pi * self.D1_av / (t_rel_rk_av * self.b_rk(self.D1_av / 2)))
        self.z_sa = int(np.pi * self.D1_av / (t_rel_sa_av * self.b_sa(self.D1_av / 2)))
        self.t_rk_av = np.pi * self.D1_av / self.z_rk
        self.t_sa_av = np.pi * self.D1_av / self.z_sa
        self.t_rel_rk_in = self.t_rk_av * self.D1_in / self.D1_av / self.b_rk(self.D1_in / 2)
        self.t_rel_sa_in = self.t_sa_av * self.D1_in / self.D1_av / self.b_sa(self.D1_in / 2)
        if self.t_rel_rk_in < self.min_t_rel_rk_in:
            self.t_rel_rk_in = self.min_t_rel_rk_in
            self.z_rk = int(np.pi * self.D1_in / (self.t_rel_rk_in * self.b_rk(self.D1_in / 2)))
            self.t_rk_av = np.pi * self.D1_av / self.z_rk
        if self.t_rel_sa_in < self.min_t_rel_sa_in:
            self.t_rel_sa_in = self.min_t_rel_sa_in
            self.z_sa = int(np.pi * self.D1_in / (self.t_rel_sa_in * self.b_sa(self.D1_in / 2)))
            self.t_sa_av = np.pi * self.D1_av / self.z_sa

    def gamma1_k_sa(self, r):
        return self.gamma1(self.alpha0(r)) * self.gamma1_k_sa_rel

    def dif_alpha0_l_gamma1_s_sa(self, r):
        dif_in = self.alpha0_l(self.D1_in / 2) - self.gamma1(self.alpha0_l(self.D1_in / 2)) * (1 - self.gamma1_k_sa_rel)
        dif_out = self.alpha0_l(self.D1_out / 2) - self.gamma1(self.alpha0_l(self.D1_out / 2)) * \
                                                   (1 - self.gamma1_k_sa_rel)
        if (dif_in < self.mindif_alpha0_l_gamma1_s_sa) and (dif_out > self.mindif_alpha0_l_gamma1_s_sa):
            dif_in = self.mindif_alpha0_l_gamma1_s_sa
            dif_out = min(self.maxdif_alpha0_l_gamma1_s_sa, dif_out)
        elif (dif_out > self.maxdif_alpha0_l_gamma1_s_sa) and (dif_in < self.maxdif_alpha0_l_gamma1_s_sa):
            dif_out = self.maxdif_alpha0_l_gamma1_s_sa
            dif_in = max(self.mindif_alpha0_l_gamma1_s_sa, dif_in)
        elif (dif_in < self.mindif_alpha0_l_gamma1_s_sa) and (dif_out < self.mindif_alpha0_l_gamma1_s_sa):
            dif_out = min(self.mindif_alpha0_l_gamma1_s_sa + abs(dif_out - dif_in), dif_out)
            dif_in = self.mindif_alpha0_l_gamma1_s_sa
        elif (dif_in > self.maxdif_alpha0_l_gamma1_s_sa) and (dif_out > self.maxdif_alpha0_l_gamma1_s_sa):
            dif_in = max(self.maxdif_alpha0_l_gamma1_s_sa - abs(dif_in - dif_out), self.mindif_alpha0_l_gamma1_s_sa)
            dif_out = self.maxdif_alpha0_l_gamma1_s_sa
        dif_int = interp1d(0.5 * np.array([self.D1_in, self.D1_out]), [dif_in, dif_out])
        return dif_int(r)

    def gamma1_s_sa(self, r):
        return self.alpha0_l(r) - self.dif_alpha0_l_gamma1_s_sa(r)

    def gamma1_k_rk(self, r):
        return self.gamma1(self.beta1(r)) * self.gamma1_k_rk_rel

    def dif_beta1_l_gamma1_s_rk(self, r):
        dif_in = self.beta1_l(self.D1_in / 2) - self.gamma1(self.beta1_l(self.D1_in / 2)) * (1 - self.gamma1_k_rk_rel)
        dif_out = self.beta1_l(self.D1_out / 2) - self.gamma1(self.beta1_l(self.D1_out / 2)) * \
                                                  (1 - self.gamma1_k_rk_rel)
        if (dif_in < self.mindif_beta1_l_gamma1_s_rk) and (dif_out > self.mindif_beta1_l_gamma1_s_rk):
            dif_in = self.mindif_beta1_l_gamma1_s_rk
            dif_out = min(self.maxdif_beta1_l_gamma1_s_rk, dif_out)
        elif (dif_out > self.maxdif_beta1_l_gamma1_s_rk) and (dif_in < self.maxdif_beta1_l_gamma1_s_rk):
            dif_out = self.maxdif_beta1_l_gamma1_s_rk
            dif_in = max(self.mindif_beta1_l_gamma1_s_rk, dif_in)
        elif (dif_in < self.mindif_beta1_l_gamma1_s_rk) and (dif_out < self.mindif_beta1_l_gamma1_s_rk):
            dif_out = min(self.mindif_beta1_l_gamma1_s_rk + abs(dif_out - dif_in), dif_out)
            dif_in = self.mindif_beta1_l_gamma1_s_rk
        elif (dif_in > self.maxdif_beta1_l_gamma1_s_rk) and (dif_out > self.maxdif_beta1_l_gamma1_s_rk):
            dif_in = max(self.maxdif_beta1_l_gamma1_s_rk - abs(dif_in - dif_out), self.mindif_beta1_l_gamma1_s_rk)
            dif_out = self.maxdif_beta1_l_gamma1_s_rk
        dif_int = interp1d(0.5 * np.array([self.D1_in, self.D1_out]), [dif_in, dif_out])
        return dif_int(r)

    def _compute_root_blade_temperature(self):
        self.T_sa_blade_in = (self.T0(0.5 * self.D1_in) + self.T1(0.5 * self.D1_in)) / 2 + \
                              0.9 * self.c0(0.5 * self.D1_in) ** 2 / (2 * self.c_p)
        self.T_rk_blade_in = (self.T1(0.5 * self.D1_in) + self.T2(0.5 * self.D1_in)) / 2 + \
                              0.9 * self.w1(0.5 * self.D1_in) ** 2 / (2 * self.c_p)

    def gamma1_s_rk(self, r):
        return self.beta1_l(r) - self.dif_beta1_l_gamma1_s_rk(r)

    def compute_sections_coordinates(self):
        self._compute_step()
        self._compute_root_blade_temperature()
        r_arr = [0.5 * (self.D1_in + i * (self.D1_out - self.D1_in)) for i in self.r_rel]
        for n, i in enumerate(r_arr):
            self._sa_sections[n].angle1 = self.alpha0_l(i)
            self._sa_sections[n].angle2 = self.alpha1_l(i)
            self._sa_sections[n].gamma1_k = self.gamma1_k_sa(i)
            self._sa_sections[n].gamma1_s = self.gamma1_s_sa(i)
            self._sa_sections[n].r = i
            self._sa_sections[n].y_s = -self._sa_sections[n].y_s
            self._sa_sections[n].y_k = -self._sa_sections[n].y_k
            self._sa_sections[n].y_av = -self._sa_sections[n].y_av
            self._rk_sections[n].angle1 = self.beta1_l(i)
            self._rk_sections[n].angle2 = self.beta2_l(i)
            self._rk_sections[n].gamma1_k = self.gamma1_k_rk(i)
            self._rk_sections[n].gamma1_s = self.gamma1_s_rk(i)
            self._rk_sections[n].r = i

    def plot2d(self, pnt_count=50, r_rel=0.5, x0=0, ymax=0.7):
        r = 0.5 * (self.D1_in + r_rel * (self.D1_out - self.D1_in))
        t_sa = self.t_sa_av * 2 * r / self.D1_av
        t_rk = self.t_rk_av * 2 * r / self.D1_av
        n_sa = int(round(ymax / t_sa)) + 1
        n_rk = int(round(ymax / t_rk)) + 1
        y0 = 0
        sa_section = BladeSection(self.alpha0_l(r), self.alpha1_l(r), self.b_a_sa, self.gamma1_s_sa(r),
                                  self.gamma1_k_sa(r), pnt_count, self.r1_sa_rel * self.b_a_sa, self.s2_sa)
        rk_section = BladeSection(self.beta1_l(r), self.beta2_l(r), self.b_a_rk, self.gamma1_s_rk(r),
                                  self.gamma1_k_rk(r), pnt_count, self.r1_rk_rel * self.b_a_rk, self.s2_rk)
        for i in range(n_sa):
            plt.plot(sa_section.x_k + x0, -sa_section.y_k + y0 + t_sa, color='red', linewidth=2)
            plt.plot(sa_section.x_s + x0, -sa_section.y_s + y0 + t_sa, color='red', linewidth=2)
            y0 += t_sa
        y0 = 0
        for i in range(n_rk):
            plt.plot(rk_section.x_k + x0 + self.b_a_sa + self.delta_a_sa, rk_section.y_k + y0, color='blue',
                     linewidth=2)
            plt.plot(rk_section.x_s + x0 + self.b_a_sa + self.delta_a_sa, rk_section.y_s + y0, color='blue',
                     linewidth=2)
            y0 += t_rk

    def plot3d(self, title, r_rel=(0, 0.5, 1), figsize=(8, 6), pnt_count=30, linewidth=2):
        r_arr = 0.5 * (self.D1_in + np.array(r_rel) * (self.D1_out - self.D1_in))
        fig = plt.figure(figsize=figsize)
        ax = Axes3D(fig)
        ax.set_title(title, fontsize=20)
        for r in r_arr:
            sa_section = BladeSection(self.alpha0_l(r), self.alpha1_l(r), self.b_a_sa, self.gamma1_s_sa(r),
                                      self.gamma1_k_sa(r), pnt_count, self.r1_sa_rel * self.b_a_sa, self.s2_sa)
            rk_section = BladeSection(self.beta1_l(r), self.beta2_l(r), self.b_a_rk, self.gamma1_s_rk(r),
                                      self.gamma1_k_rk(r), pnt_count, self.r1_rk_rel * self.b_a_rk, self.s2_rk)
            ax.plot(xs=sa_section.x_s, ys=-sa_section.y_s + self.t_sa_av, zs=r, color='red', linewidth=linewidth)
            ax.plot(xs=sa_section.x_k, ys=-sa_section.y_k + self.t_sa_av, zs=r, color='red', linewidth=linewidth)
            ax.plot(xs=rk_section.x_s + self.b_a_sa + self.delta_a_sa, ys=rk_section.y_s, zs=r, color='blue',
                    linewidth=linewidth)
            ax.plot(xs=rk_section.x_k + self.b_a_sa + self.delta_a_sa, ys=rk_section.y_k, zs=r, color='blue',
                    linewidth=linewidth)
        plt.show()


def get_stage_profiling_object(stage_geom: StageGeomAndHeatDrop, stage_gas_dynamics: StageGasDynamics,
                               prof_type: ProfilingType, p0_stag, c0, alpha0, pnt_count=20, r_rel=(0, 0.5, 1),
                               gamma1_k_sa_rel=0.25, gamma1_k_rk_rel=0.25, s2_sa=0.001, s2_rk=0.001) -> StageProfiling:
    """

    :param stage_geom:
    :param stage_gas_dynamics:
    :param prof_type:
    :param alpha0: угол входа в СА в зависимости от радиуса
    :param pnt_count:  количесиство точек для нахождения координат сторон профиля в каждом сечении
    :param r_rel: относительные радиусы сечений для нахождения координат
    :param gamma1_k_sa_rel: отношение углов gamma1_k и gamma1 для СА
    :param gamma1_k_rk_rel: отношение углов gamma1_k и gamma1 для РК
    :param s2_sa: толщина выходной кромки СА
    :param s2_rk: толщина выходной кромки РК
    :return: StageProfiling
    """
    result = StageProfiling(prof_type, p0_stag, stage_gas_dynamics.T0_stag, c0,
                            stage_gas_dynamics.c_p_gas, stage_gas_dynamics.k_gas,
                            stage_geom.D1 - stage_geom.l1, stage_geom.D1, stage_geom.D1 + stage_geom.l1,
                            stage_gas_dynamics.n, stage_gas_dynamics.c1, stage_gas_dynamics.alpha1,
                            stage_gas_dynamics.T_st_stag, stage_gas_dynamics.L_u, stage_gas_dynamics.c2_a,
                            stage_gas_dynamics.c2_u, alpha0, pnt_count, stage_geom.b_sa, stage_geom.b_rk,
                            stage_geom.delta_a_sa, r_rel=r_rel, gamma1_k_sa_rel=gamma1_k_sa_rel,
                            gamma1_k_rk_rel=gamma1_k_rk_rel, s2_sa=s2_sa, s2_rk=s2_rk)
    return result


class TurbineProfiling:
    def __init__(self, turbine: Turbine, profiling_type: ProfilingType, pnt_count=20, r_rel=(0, 0.5, 1)):
        self._turbine = turbine
        self._stages_profile = []
        self._profiling_type = profiling_type
        self._pnt_count = pnt_count
        self._r_rel = r_rel
        for i in range(len(turbine)):
            deg = np.pi / 180
            if i == 0:
                stage_pr = get_stage_profiling_object(turbine.geom[i], turbine[i], profiling_type,
                                                      lambda r: turbine[i].p0_stag, lambda r: 100, lambda r: 90 * deg,
                                                      pnt_count=pnt_count, r_rel=r_rel)
                self._stages_profile.append(stage_pr)
            else:
                stage_pr = get_stage_profiling_object(turbine.geom[i], turbine[i], profiling_type, self[i - 1].p2_stag,
                                                      self[i - 1].c0, self[i - 1].alpha2, pnt_count=pnt_count,
                                                      r_rel=r_rel)
                self._stages_profile.append(stage_pr)

    @property
    def r_rel(self) -> tuple:
        return self._r_rel

    @r_rel.setter
    def r_rel(self, value: tuple):
        self._r_rel = value
        for i in self:
            i.r_rel = value

    @property
    def pnt_count(self):
        return self._pnt_count

    @pnt_count.setter
    def pnt_count(self, value):
        self._pnt_count = value
        for i in self:
            i.pnt_count = value

    @property
    def profiling_type(self) -> ProfilingType:
        return self._profiling_type

    @profiling_type.setter
    def profiling_type(self, value: ProfilingType):
        self._profiling_type = value
        for i in self:
            i.profiling_type = value

    def __getitem__(self, item) -> StageProfiling:
        if 0 <= item < self._turbine.stage_number:
            return self._stages_profile[item]
        else:
            raise IndexError('invalid index')

    def __len__(self):
        return self._turbine.stage_number

    def __iter__(self):
        self._num = 0
        return self

    def __next__(self) -> StageProfiling:
        if self._num < self._turbine.stage_number:
            current = self._stages_profile[self._num]
            self._num += 1
            return current
        else:
            raise StopIteration()

    def set_gamma1_k_rel(self, gamma1_k_sa_rel, gamma1_k_rk_rel):
        """Задает величину gamma1_k / gamma1 для СА и РК всех ступеней"""
        for i in self:
            i.gamma1_k_sa_rel = gamma1_k_sa_rel
            i.gamma1_k_rk_rel = gamma1_k_rk_rel

    def set_edge_thickness(self, s2_sa, s2_rk):
        """Задает толщину выходных кромок лопаток СА и РК на всех ступенях"""
        for i in self:
            i.s2_sa = s2_sa
            i.s2_rk = s2_rk

    def set_inlet_edge_rel_radius(self, r1_rk_rel, r1_sa_rel):
        """Задает относительный радиус закругления входных кромок лопаток СА и РК на всех ступенях"""
        for i in self:
            i.r1_rk_rel = r1_rk_rel
            i.r1_sa_rel = r1_sa_rel

    def compute_profile(self):
        for i in self:
            i.compute_sections_coordinates()

    def save(self):
        file = open(os.path.join(os.path.dirname(__file__), 'profiling_results'), 'wb')
        arr = [dict() for _ in range(len(self))]
        for n1, i1 in enumerate(self):
            arr[n1]['sa'] = dict()
            arr[n1]['sa']['z'] = i1.z_sa
            arr[n1]['sa']['sections'] = []
            arr[n1]['sa']['D0_in'] = self._turbine.geom[n1].D0 - self._turbine.geom[n1].l0
            arr[n1]['sa']['D05_in'] = self._turbine.geom[n1].D05 - self._turbine.geom[n1].l05
            arr[n1]['sa']['D0_out'] = self._turbine.geom[n1].D0 + self._turbine.geom[n1].l0
            arr[n1]['sa']['D05_out'] = self._turbine.geom[n1].D05 + self._turbine.geom[n1].l05
            arr[n1]['sa']['b_a'] = i1.b_a_sa
            arr[n1]['rk'] = dict()
            arr[n1]['rk']['z'] = i1.z_rk
            arr[n1]['rk']['sections'] = []
            arr[n1]['rk']['D1_in'] = self._turbine.geom[n1].D1 - self._turbine.geom[n1].l1
            arr[n1]['rk']['D2_in'] = self._turbine.geom[n1].D2 - self._turbine.geom[n1].l2
            arr[n1]['rk']['D1_out'] = self._turbine.geom[n1].D1 + self._turbine.geom[n1].l1
            arr[n1]['rk']['D2_out'] = self._turbine.geom[n1].D2 + self._turbine.geom[n1].l2
            arr[n1]['rk']['delta_r'] = self._turbine.geom[n1].delta_r_rk
            arr[n1]['rk']['b_a'] = i1.b_a_rk
            arr[n1]['rk']['delta_a_sa'] = self._turbine.geom[n1].delta_a_sa
            arr[n1]['rk']['delta_a_rk'] = self._turbine.geom[n1].delta_a_rk
            arr[n1]['rk']['T_rk_blade_in'] = i1.T_rk_blade_in
            for i2 in i1:
                arr[n1]['sa']['sections'].append(i2['sa'])
                arr[n1]['rk']['sections'].append(i2['rk'])
        pk.dump(arr, file)
        file.close()

    def plot2d(self, figsize=(10, 8), pnt_count=50, r_rel=0.5, ymax=0.7):
        plt.figure(figsize=figsize)
        plt.title(r'$r_{rel} = %s$' % r_rel, fontsize=20)
        x0 = 0
        for i in range(len(self)):
            self[i].plot2d(pnt_count, r_rel, x0, ymax)
            x0 += self._turbine.geom[i].b_sa + self._turbine.geom[i].b_rk + self._turbine.geom[i].delta_a_sa + \
                  self._turbine.geom[i].delta_a_rk
        plt.ylim(0, ymax)
        plt.grid()
        plt.show()


if __name__ == '__main__':
    deg = np.pi / 180
    turbine = Turbine(TurbineType.Compressor, gamma_av=4 * deg, gamma_sum=10 * deg)
    turbine.alpha11 = 18 * deg
    turbine.alpha_air = 2.87
    turbine.c21_init = 220
    turbine.eta_t_stag_cycle = 0.91
    turbine.G_turbine = 30
    turbine.H01_init = 60e3
    turbine.H_t_stag_cycle = 150e3
    turbine.k_n = 6.8
    turbine.l1_D1_ratio = 1 / 3
    turbine.L_t_cycle = 130e3
    turbine.n = 10e3
    turbine.p_g_stag = 11e5
    turbine.T_g_stag = 1400
    turbine.T_t_stag_cycle = 950
    turbine.p_t_stag_cycle = 850e3
    turbine.phi1 = 0.97
    turbine.rho1 = 0.3
    turbine.stage_number = 2
    turbine.set_delta_a_b_ratio(0.22, 0)
    turbine.set_l_b_ratio(1.8, 0.2, 0.9)
    turbine.compute_geometry(compute_heat_drop_auto=False)
    turbine.geom[0].H0 = 60e3
    turbine.geom[1].H0 = 200e3
    turbine.geom[0].rho = 0.45
    turbine.geom[1].rho = 0.35
    turbine.compute_stages_gas_dynamics()
    turbine[0].plot_velocity_triangle()
    turbine[1].plot_velocity_triangle()
    print(turbine.geom.sigma_l / 1e6)
    # print(turbine[0].alpha2 / deg)
    # print(turbine[0].H0)
    # print(turbine[1].alpha2 / deg)
    # print(turbine[1].H0)
    turbine_profiling = TurbineProfiling(turbine, ProfilingType.ConstantAngle)
    turbine_profiling.set_gamma1_k_rel(0.25, 0.35)
    turbine_profiling.compute_profile()
    turbine_profiling.plot2d(r_rel=0.5)
    # for i in turbine_profiling:
    #     print(i.z_sa, i.b_sa, i.gamma_sa / deg, i.t_rel_sa, i.gamma1_s_sa(0.5 * i.D1_in) / deg,
    #           i.gamma1_k_sa(0.5 * i.D1_in) / deg)
    #     print(i.z_rk, i.b_rk, i.gamma_rk / deg, i.t_rel_rk, i.gamma1_s_rk(i.D1_av / 2) / deg,
    #           i.beta1_l(0.5 * i.D1_av) / deg)
    turbine_profiling[0].plot_velocity_triangles()
    # turbine_profiling[0].plot3d(title='Stage 1', r_rel=np.linspace(0, 1, 20))
    # turbine_profiling[1].plot3d(title='Stage 2', r_rel=np.linspace(0, 1, 20))
