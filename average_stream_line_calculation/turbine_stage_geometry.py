from gas import KeroseneCombustionProducts
import functions as func
import logging
import matplotlib.pyplot as plt
from gas_dynamics import *
from scipy.optimize import fsolve
import os


log_filename = os.path.join(os.path.dirname(__file__), 'average_streamline_calculation.log')
logger = func.create_logger(__name__, logging.DEBUG, filename=log_filename, filemode='a',
                            add_console_handler=False)


class InvalidStageSizeValue(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)
        self.message = message


class StageGeomAndHeatDrop:
    def __init__(self):
        self._gamma_out = None
        self._gamma_in = None
        self._gamma_av = None
        self._l1_b_sa_ratio = None
        self._l2_b_rk_ratio = None
        self._delta_a_b_sa_ratio = 0.22
        self._delta_a_b_rk_ratio = 0.22
        self._mu = 1
        self._k_b_sa = None
        self._k_b_rk = None
        self._k_delta_a_sa = None
        self._k_delta_a_rk = None
        self._l1 = None
        self._D1 = None
        self._delta_r_rk_l2_ratio = 0.01
        self._n = None
        self._x0 = None
        self._rho = None
        self._phi = 0.97
        self._psi = 0.97
        self._epsilon = 1
        self._g_lk = 0
        self._g_ld = 0
        self._g_lb = 0
        self.D2 = None
        self.D0 = None
        self.D05 = None
        self.l0 = None
        self.l05 = None
        self.l2 = None
        self.u1 = None
        self.u2 = None
        self.u_av = None
        self.H0 = None
        self.delta_a_sa = None
        self.delta_a_rk = None
        self.b_sa = None
        self.b_rk = None
        self.length = None

    def str(self):
        str_arr = str(self).split()
        return str_arr[0][1:] + ' ' + str_arr[1]

    def _compute_coefficients(self):
        logger.info('%s _compute_coefficients' % (self.str()))
        self._k_b_sa = 1 / self._l1_b_sa_ratio
        self._k_b_rk = 1 / self._l2_b_rk_ratio
        self._k_delta_a_sa = self._delta_a_b_sa_ratio / self._l1_b_sa_ratio
        self._k_delta_a_rk = self._delta_a_b_rk_ratio / self._l2_b_rk_ratio

    def _compute_geometry(self):
        logger.info('%s _compute_geometry' % (self.str()))
        self.b_sa = self.l1 * self.k_b_sa
        self.delta_a_sa = self.l1 * self.k_delta_a_sa
        self.l0 = self.l1 - (np.tan(self.gamma_in) + np.tan(self.gamma_out)) * (self.b_sa + self.delta_a_sa)
        self.l2 = self.l1 / (1 - self.k_b_rk * (np.tan(self.gamma_in) + np.tan(self.gamma_out)))
        self.delta_r_rk = self.delta_r_rk_l2_ratio * self.l2
        self.b_rk = self.l2 * self.k_b_rk
        self.delta_a_rk = self.l2 * self.k_delta_a_rk
        self.D2 = self.D1 + 2 * np.tan(self.gamma_av) * self.b_rk
        self.D05 = self.D1 - 2 * np.tan(self.gamma_av) * self.delta_a_sa
        self.D0 = self.D05 - 2 * np.tan(self.gamma_av) * self.b_sa
        self.l05 = self.l0 + self.b_sa * (np.tan(self.gamma_in) + np.tan(self.gamma_out))
        self.length = self.b_sa + self.delta_a_sa + self.b_rk + self.delta_a_rk
        self.l0_next = self.l2 + self.delta_a_rk * (np.tan(self.gamma_in) + np.tan(self.gamma_out))
        self.A1 = np.pi * self.D1 * self.l1
        self.A2 = np.pi * self.D2 * self.l2

    def _compute_velocities(self):
        logger.info('%s _compute_velocities' % (self.str()))
        self.u1 = np.pi * self.D1 * self.n / 60
        self.u2 = np.pi * self.D2 * self.n / 60
        self.u_av = np.pi * 0.5 * (self.D2 + self.D1) * self.n / 60

    def _computing_geometry_condition(self):
        return not ((self._l2_b_rk_ratio is None) or (self._l1_b_sa_ratio is None) or
                    (self._delta_a_b_rk_ratio is None) or (self._delta_a_b_sa_ratio is None) or
                    (self._gamma_av is None) or (self._gamma_out is None) or (self._gamma_in is None) or
                    (self._l1 is None) or (self._D1 is None) or (self._delta_r_rk_l2_ratio is None))

    def _computing_coefficients_condition(self):
        return not ((self._l2_b_rk_ratio is None) or (self._l1_b_sa_ratio is None) or
                    (self._delta_a_b_rk_ratio is None) or (self._delta_a_b_sa_ratio is None))

    def _computing_velocities_condition(self):
        return not ((self._l2_b_rk_ratio is None) or (self._l1_b_sa_ratio is None) or
                    (self._delta_a_b_rk_ratio is None) or (self._delta_a_b_sa_ratio is None) or
                    (self._gamma_av is None) or (self._gamma_out is None) or (self._gamma_in is None) or
                    (self._l1 is None) or (self._D1 is None) or (self._delta_r_rk_l2_ratio is None) or
                    (self._n is None))

    def plot(self):
        logger.info('%s plot' % (self.str()))
        x_sa_arr = np.array([self.x0, self.x0, self.x0 + self.b_sa, self.x0 + self.b_sa, self.x0])
        y_sa_arr = np.array([0.5 * (self.D0 - self.l0), 0.5 * (self.D0 + self.l0), 0.5 * (self.D05 + self.l05),
                            0.5 * (self.D05 - self.l05), 0.5 * (self.D0 - self.l0)])
        x0_rk = self.x0 + self.b_sa + self.delta_a_sa
        x_rk_arr = np.array([x0_rk, x0_rk, x0_rk + self.b_rk, x0_rk + self.b_rk, x0_rk])
        y_rk_arr = np.array([0.5 * (self.D1 - self.l1), 0.5 * (self.D1 + self.l1) - self.delta_r_rk,
                             0.5 * (self.D2 + self.l2) - self.delta_r_rk, 0.5 * (self.D2 - self.l2),
                             0.5 * (self.D1 - self.l1)])
        x_out_arr = np.array([self.x0, x0_rk + self.b_rk + self.delta_a_rk])
        y_out_arr = np.array([0.5 * (self.D0 + self.l0),
                              0.5 * (self.D0 + self.l0) + np.tan(self.gamma_out) * self.length])
        y_av_arr = np.array([0.5 * self.D0, 0.5 * self.D0 + np.tan(self.gamma_av) * self.length])
        plt.plot(x_sa_arr, y_sa_arr, linewidth=2, color='black')
        plt.plot(x_rk_arr, y_rk_arr, linewidth=2, color='black')
        plt.plot(x_out_arr, y_out_arr, linewidth=2, color='black')
        plt.plot(x_out_arr, y_av_arr, '--', linewidth=2, color='black')
        plt.plot(x_sa_arr, -y_sa_arr, linewidth=2, color='black')
        plt.plot(x_rk_arr, -y_rk_arr, linewidth=2, color='black')
        plt.plot(x_out_arr, -y_out_arr, linewidth=2, color='black')
        plt.plot(x_out_arr, -y_av_arr, '--', linewidth=2, color='black')

    @property
    def g_ld(self):
        assert self._g_ld is not None, 'g_ld must not be None'
        return self._g_ld

    @g_ld.setter
    def g_ld(self, value):
        self._g_ld = value

    @property
    def g_lb(self):
        assert self._g_lb is not None, 'g_lb must not be None'
        return self._g_lb

    @g_lb.setter
    def g_lb(self, value):
        self._g_lb = value

    @property
    def g_lk(self):
        assert self._g_lk is not None, 'g_lk must not be None'
        return self._g_lk

    @g_lk.setter
    def g_lk(self, value):
        self._g_lk = value

    @property
    def epsilon(self):
        assert self._epsilon is not None, 'epsilon must not be None'
        return self._epsilon

    @epsilon.setter
    def epsilon(self, value):
        self._epsilon = value

    @property
    def psi(self):
        assert self._psi is not None, 'psi must not be None'
        return self._psi

    @psi.setter
    def psi(self, value):
        self._psi = value

    @property
    def phi(self):
        assert self._phi is not None, 'phi must not be None'
        return self._phi

    @phi.setter
    def phi(self, value):
        self._phi = value

    @property
    def rho(self):
        assert self._rho is not None, 'rho must not be None'
        return self._rho

    @rho.setter
    def rho(self, value):
        self._rho = value

    @property
    def x0(self):
        assert self._x0 is not None, 'x0 must not be None'
        return self._x0

    @x0.setter
    def x0(self, value):
        self._x0 = value

    @property
    def delta_r_rk_l2_ratio(self):
        assert self._delta_r_rk_l2_ratio is not None, 'delta_r_rk_l2_ratio must not be None'
        return self._delta_r_rk_l2_ratio

    @delta_r_rk_l2_ratio.setter
    def delta_r_rk_l2_ratio(self, value):
        self._delta_r_rk_l2_ratio = value

    @property
    def n(self):
        assert self._n is not None, 'n_rel must not be None'
        return self._n

    @n.setter
    def n(self, value):
        self._n = value
        if self._computing_velocities_condition():
            self._compute_velocities()

    @property
    def gamma_out(self):
        assert self._gamma_out is not None, 'gamma_out must not be None'
        return self._gamma_out

    @gamma_out.setter
    def gamma_out(self, value):
        self._gamma_out = value
        if self._computing_geometry_condition():
            self._compute_geometry()
        if self._computing_velocities_condition():
            self._compute_velocities()

    @property
    def gamma_in(self):
        assert self._gamma_in is not None, 'gamma_in must not be None'
        return self._gamma_in

    @gamma_in.setter
    def gamma_in(self, value):
        self._gamma_in = value
        if self._computing_geometry_condition():
            self._compute_geometry()
        if self._computing_velocities_condition():
            self._compute_velocities()

    @property
    def gamma_av(self):
        assert self._gamma_av is not None, 'gamma_av must not be None'
        return self._gamma_av

    @gamma_av.setter
    def gamma_av(self, value):
        self._gamma_av = value
        if self._computing_geometry_condition():
            self._compute_geometry()
        if self._computing_velocities_condition():
            self._compute_velocities()

    @property
    def D1(self):
        assert self._D1 is not None, 'D1 must not be None'
        return self._D1

    @D1.setter
    def D1(self, value):
        self._D1 = value
        if self._computing_geometry_condition():
            self._compute_geometry()
        if self._computing_velocities_condition():
            self._compute_velocities()

    @property
    def l1(self):
        assert self._l1 is not None, 'l1 must not be None'
        return self._l1

    @l1.setter
    def l1(self, value):
        self._l1 = value
        if self._computing_geometry_condition():
            self._compute_geometry()
        if self._computing_velocities_condition():
            self._compute_velocities()

    @property
    def k_delta_a_rk(self):
        return self._k_delta_a_rk

    @property
    def k_delta_a_sa(self):
        return self._k_delta_a_sa

    @property
    def k_b_rk(self):
        return self._k_b_rk

    @property
    def k_b_sa(self):
        return self._k_b_sa

    @property
    def mu(self):
        assert self._mu is not None, 'mu must not be None'
        return self._mu

    @mu.setter
    def mu(self, value):
        self._mu = value

    @property
    def l1_b_sa_ratio(self):
        assert self._l1_b_sa_ratio is not None, 'l_b_sa_ratio must not be None'
        return self._l1_b_sa_ratio

    @l1_b_sa_ratio.setter
    def l1_b_sa_ratio(self, value):
        self._l1_b_sa_ratio = value
        if self._computing_coefficients_condition():
            self._compute_coefficients()
        if self._computing_geometry_condition():
            self._compute_geometry()
        if self._computing_velocities_condition():
            self._compute_velocities()

    @property
    def l2_b_rk_ratio(self):
        assert self._l2_b_rk_ratio is not None, 'l_b_rk_ratio must not be None'
        return self._l2_b_rk_ratio

    @l2_b_rk_ratio.setter
    def l2_b_rk_ratio(self, value):
        self._l2_b_rk_ratio = value
        if self._computing_coefficients_condition():
            self._compute_coefficients()
        if self._computing_geometry_condition():
            self._compute_geometry()
        if self._computing_velocities_condition():
            self._compute_velocities()

    @property
    def delta_a_b_sa_ratio(self):
        assert self._delta_a_b_sa_ratio is not None, 'delta_a_b_sa_ratio must not be None'
        return self._delta_a_b_sa_ratio

    @delta_a_b_sa_ratio.setter
    def delta_a_b_sa_ratio(self, value):
        self._delta_a_b_sa_ratio = value
        if self._computing_coefficients_condition():
            self._compute_coefficients()
        if self._computing_geometry_condition():
            self._compute_geometry()
        if self._computing_velocities_condition():
            self._compute_velocities()

    @property
    def delta_a_b_rk_ratio(self):
        assert self._delta_a_b_rk_ratio is not None, 'delta_a_rk must not be None'
        return self._delta_a_b_rk_ratio

    @delta_a_b_rk_ratio.setter
    def delta_a_b_rk_ratio(self, value):
        self._delta_a_b_rk_ratio = value
        if self._computing_coefficients_condition():
            self._compute_coefficients()
        if self._computing_geometry_condition():
            self._compute_geometry()
        if self._computing_velocities_condition():
            self._compute_velocities()


class TurbineGeomAndHeatDropDistribution:
    def __init__(self, stage_number, eta_t_stag, H_t_stag, c21, n_rel, work_fluid: KeroseneCombustionProducts, T_g_stag,
                 p_g_stag, alpha_air, G_turbine, l1_D1_ratio, H01, rho1, phi1, alpha11, k_n, sigma_l, T_t_stag,
                 p_t_stag, **kwargs):
        """
        :param stage_number:
        :param eta_t_stag:
        :param H_t_stag:
        :param c21:
        :param n_rel: отношение частоты вращения к максимально возможной
        :param work_fluid:
        :param T_g_stag:
        :param p_g_stag:
        :param alpha_air: коэффициент избытка воздуха
        :param G_turbine:
        :param l1_D1_ratio:
        :param H01:
        :param rho1:
        :param phi1:
        :param alpha11: угол потока после СА первой ступени
        :param k_n:
        :param sigma_l:
        :param T_t_stag:
        :param p_t_stag:
        :param kwargs: gamma_av, gamma_sum, gamma_in, gamma_out
        """
        self._stage_number = stage_number
        self._eta_t_stag = eta_t_stag
        self._H_t_stag = H_t_stag
        self._c21 = c21
        self._n_rel = n_rel
        self._work_fluid = work_fluid
        self._T_g_stag = T_g_stag
        self._p_g_stag = p_g_stag
        self._alpha_air = alpha_air
        self._G_turbine = G_turbine
        self._l1_D1_ratio = l1_D1_ratio
        self._H01 = H01
        self._rho1 = rho1
        self._phi1 = phi1
        self._alpha11 = alpha11
        self._k_n = k_n
        self._sigma_l = sigma_l
        self._T_t_stag = T_t_stag
        self._p_t_stag = p_t_stag
        self._kwargs = kwargs
        if ('gamma_av' in kwargs) and ('gamma_sum' in kwargs):
            self._gamma_av = kwargs['gamma_av']
            self._gamma_in = None
            self._gamma_out = None
            self._gamma_sum = kwargs['gamma_sum']
        elif ('gamma_in' in kwargs) and ('gamma_out' in kwargs):
            self._gamma_av = None
            self._gamma_in = kwargs['gamma_in']
            self._gamma_out = kwargs['gamma_out']
            self._gamma_sum = None
        else:
            assert False, 'gamma_av and gamma_sum or gamma_in and gamma_out must be set'
        self._stages = [StageGeomAndHeatDrop() for _ in range(self._stage_number)]

    @property
    def p_t_stag(self):
        assert self._p_t_stag is not None, 'p_t_stag must not be None'
        return self._p_t_stag

    @p_t_stag.setter
    def p_t_stag(self, value):
        self._p_t_stag = value

    @property
    def T_t_stag(self):
        assert self._T_t_stag is not None, 'T_t_stag must not be None'
        return self._T_t_stag

    @T_t_stag.setter
    def T_t_stag(self, value):
        self._T_t_stag = value

    @property
    def sigma_l(self):
        assert self._sigma_l is not None, 'sigma_l must not be None'
        return self._sigma_l

    @sigma_l.setter
    def sigma_l(self, value):
        self._sigma_l = value

    @property
    def k_n(self):
        assert self._k_n is not None, 'k_n must not be None'
        return self._k_n

    @k_n.setter
    def k_n(self, value):
        self._k_n = value

    @property
    def gamma_sum(self):
        if 'gamma_sum' in self._kwargs:
            assert self._gamma_sum is not None, 'gamma_sum must not be None'
        return self._gamma_sum

    @gamma_sum.setter
    def gamma_sum(self, value):
        self._gamma_sum = value

    @property
    def gamma_out(self):
        if 'gamma_out' in self._kwargs:
            assert self._gamma_out is not None, 'gamma_out must not be None'
        return self._gamma_out

    @gamma_out.setter
    def gamma_out(self, value):
        self._gamma_out = value

    @property
    def gamma_in(self):
        if 'gamma_in' in self._kwargs:
            assert self._gamma_in is not None, 'gamma_in must not be None'
        return self._gamma_in

    @gamma_in.setter
    def gamma_in(self, value):
        self._gamma_in = value

    @property
    def gamma_av(self):
        if 'gamma_av' in self._kwargs:
            assert self._gamma_av is not None, 'gamma_av must not be None'
        return self._gamma_av

    @gamma_av.setter
    def gamma_av(self, value):
        self._gamma_av = value

    @property
    def alpha11(self):
        assert self._alpha11 is not None, 'alpha11 must not be None'
        return self._alpha11

    @alpha11.setter
    def alpha11(self, value):
        self._alpha11 = value

    @property
    def phi1(self):
        assert self._phi1 is not None, 'phi1 must not be None'
        return self._phi1

    @phi1.setter
    def phi1(self, value):
        self._phi1 = value

    @property
    def rho1(self):
        assert self._rho1 is not None, 'rho1 must not be None'
        return self._rho1

    @rho1.setter
    def rho1(self, value):
        self._rho1 = value

    @property
    def H01(self):
        assert self._H01 is not None, 'H01 must not be None'
        return self._H01

    @H01.setter
    def H01(self, value):
        self._H01 = value

    @property
    def l1_D1_ratio(self):
        assert self._l1_D1_ratio is not None, 'l1_D1_ratio must not be None'
        return self._l1_D1_ratio

    @l1_D1_ratio.setter
    def l1_D1_ratio(self, value):
        self._l1_D1_ratio = value

    @property
    def G_turbine(self):
        assert self._G_turbine is not None, 'G_stage_in must not be None'
        return self._G_turbine

    @G_turbine.setter
    def G_turbine(self, value):
        self._G_turbine = value

    @property
    def alpha_air(self):
        assert self._alpha_air is not None, 'alpha_air must not be None'
        return self._alpha_air

    @alpha_air.setter
    def alpha_air(self, value):
        self._alpha_air = value

    @property
    def p_g_stag(self):
        assert self._p_g_stag is not None, 'p_g_stag must not be None'
        return self._p_g_stag

    @p_g_stag.setter
    def p_g_stag(self, value):
        self._p_g_stag = value

    @property
    def T_g_stag(self):
        assert self._T_g_stag is not None, 'T_g_stag must not be None'
        return self._T_g_stag

    @T_g_stag.setter
    def T_g_stag(self, value):
        self._T_g_stag = value

    @property
    def work_fluid(self) -> KeroseneCombustionProducts:
        assert self._work_fluid is not None, 'work_fluid must not be None'
        return self._work_fluid

    @work_fluid.setter
    def work_fluid(self, value: KeroseneCombustionProducts):
        self._work_fluid = value

    @property
    def n_rel(self):
        assert self._n_rel is not None, 'n_rel must not be None'
        return self._n_rel

    @n_rel.setter
    def n_rel(self, value):
        self._n_rel = value

    @property
    def c21(self):
        assert self._c21 is not None, 'c21 must not be None'
        return self._c21

    @c21.setter
    def c21(self, value):
        self._c21 = value

    @property
    def H_t_stag(self):
        assert self._H_t_stag is not None, 'H_t_stag must not be None'
        return self._H_t_stag

    @H_t_stag.setter
    def H_t_stag(self, value):
        self._H_t_stag = value

    @property
    def eta_t_stag(self):
        assert self._eta_t_stag is not None, 'eta_t_stag must not be None'
        return self._eta_t_stag

    @eta_t_stag.setter
    def eta_t_stag(self, value):
        self._eta_t_stag = value

    @property
    def stage_number(self):
        assert self._stage_number is not None, 'stage_number must not be None'
        return self._stage_number

    @stage_number.setter
    def stage_number(self, value):
        self._stage_number = value
        self._stages = [StageGeomAndHeatDrop() for _ in range(self._stage_number)]

    @property
    def last(self) -> StageGeomAndHeatDrop:
        return self._stages[self.stage_number - 1]

    @property
    def first(self) -> StageGeomAndHeatDrop:
        return self._stages[0]

    def str(self):
        str_arr = str(self).split()
        return str_arr[0][1:] + ' ' + str_arr[1]

    def _compute_d1_and_l1(self):
        logger.info('%s _compute_d1_and_l1' % (self.str()))
        self.work_fluid.__init__()
        self.dT11_rel = 1
        self._iter_number_d1_l1 = 0
        while self.dT11_rel >= 0.001:
            self._iter_number_d1_l1 += 1
            self._compute_t_11()
        self.p_11 = self.p_g_stag * (1 - self.H_s1 /
                                     (self.c_p_gas11 * self.T_g_stag)) ** (self.k_gas11 / (self.k_gas11 - 1))
        self.rho11 = self.p_11 / (self.work_fluid.R * self.T_11)
        self.D1_l1_product = self.G_turbine / (np.pi * self.rho11 * np.sin(self.alpha11) * self.c11)
        logger.debug('%s _compute_d1_and_l1 D1_l1_product = %s' % (self.str(), self.D1_l1_product))
        self.D1 = np.sqrt(self.D1_l1_product / self.l1_D1_ratio)
        self.l1 = self.D1_l1_product / self.D1
        logger.debug('%s _compute_d1_and_l1 D1 = %s, l1 = %s' % (self.str(), self.D1, self.l1))

    def _compute_t_11(self):
        """Вычисляет величину T_11"""
        logger.info('%s _compute_t_11' % (self.str()))
        logger.debug('%s _compute_t_11 iter_number = %s' % (self.str(), self._iter_number_d1_l1))
        self.work_fluid.T1 = self.T_g_stag
        self.work_fluid.alpha = self.alpha_air
        self.H_s1 = self.H01 * (1 - self.rho1)
        self.c11 = self.phi1 * np.sqrt(2 * self.H_s1)
        self.k_gas11 = self.work_fluid.k_av_int
        self.c_p_gas11 = self.work_fluid.c_p_av_int
        logger.debug('%s _compute_t_11 c_p_gas11 = %s' % (self.str(), self.k_gas11))
        self.T_11 = self.T_g_stag - self.H_s1 * self.phi1 ** 2 / self.c_p_gas11
        self.dT11_rel = abs(self.T_11 - self.work_fluid.T2) / self.work_fluid.T2
        self.work_fluid.T2 = self.T_11
        logger.debug('%s _compute_t_11 dT11_rel = %s' % (self.str(), self.dT11_rel))

    def _compute_angles(self):
        logger.info('%s _compute_angles' % (self.str()))
        if ('gamma_in' in self._kwargs) and ('gamma_out' in self._kwargs):
            self.gamma_av = np.arctan(0.5 * (np.tan(self.gamma_out) - np.tan(self.gamma_in)))
            self.gamma_sum = self.gamma_in + self.gamma_out
        elif ('gamma_av' in self._kwargs) and ('gamma_sum' in self._kwargs):
            a = np.tan(self.gamma_sum)
            b = 2 - 2 * np.tan(self.gamma_sum) * np.tan(self.gamma_av)
            c = -2 * np.tan(self.gamma_av) - np.tan(self.gamma_sum)
            d = b**2 - 4 * a * c
            if d < 0:
                raise InvalidStageSizeValue('d < 0')
            self.gamma_out = np.arctan((-b + np.sqrt(d)) / (2 * a))
            self.gamma_in = np.arctan(np.tan(self.gamma_out) - 2 * np.tan(self.gamma_av))

    def __getitem__(self, item) -> StageGeomAndHeatDrop:
        if 0 <= item < self._stage_number:
            return self._stages[item]
        else:
            raise IndexError('invalid index')

    def __len__(self):
        return len(self._stages)

    def __iter__(self):
        self._num = 0
        return self

    def __next__(self) -> StageGeomAndHeatDrop:
        if self._num < self.stage_number:
            current = self._stages[self._num]
            self._num += 1
            return current
        else:
            raise StopIteration()

    def _compute_linear_dimensions(self):
        logger.info('%s _compute_linear_dimensions' % (self.str()))
        for num, item in enumerate(self._stages):
            item.gamma_av = self.gamma_av
            item.gamma_out = self.gamma_out
            item.gamma_in = self.gamma_in
            if num == 0:
                item.l1 = self.l1
                item.D1 = self.D1
            else:
                l0 = self._stages[num - 1].l0_next
                item.l1 = l0 / (1 - (item.k_b_sa + item.k_delta_a_sa) *
                                (np.tan(self.gamma_in) + np.tan(self.gamma_out)))
                b_sa = item.l1 * item.k_b_sa
                delta_a_sa = item.l1 * item.k_delta_a_sa
                item.D1 = self._stages[num - 1].D2 + 2 * np.tan(self.gamma_av) * (self._stages[num - 1].delta_a_rk +
                                                                                  b_sa + delta_a_sa)

    def _compute_geometry(self):
        logger.info('%s  _compute_geometry' % (self.str()))
        self._compute_d1_and_l1()
        self._compute_angles()
        self._compute_linear_dimensions()

    def plot_geometry(self, figsize=(5, 11), title='Turbine geometry'):
        logger.info('%s  plot_geometry' % (self.str()))
        plt.figure(figsize=figsize)
        for num, item in enumerate(self._stages):
            if num == 0:
                item.x0 = 0
                item.plot()
            else:
                item.x0 = self._stages[num - 1].x0 + self._stages[num - 1].length
                item.plot()
        plt.grid()
        plt.title(title, fontsize=20)
        plt.xlim(-0.01, self._stages[self.stage_number - 1].x0 + self._stages[self.stage_number - 1].length + 0.01)
        plt.show()

    def _compute_heat_drop_distribution(self):
        logger.info('%s  _compute_heat_drop_distribution' % (self.str()))
        u_av_squared_sum = 0
        for num, item in enumerate(self._stages):
            u_av_squared_sum += item.u_av ** 2
        for num, item in enumerate(self._stages):
            if num == 0:
                item.H0 = self.H_t * (1 + self.alpha) * item.u_av**2 / u_av_squared_sum + 0.5 * (item.mu * self.c21)**2
            else:
                item.H0 = self.H_t * (1 + self.alpha) * item.u_av ** 2 / u_av_squared_sum

    def plot_heat_drop_distribution(self, figsize=(9, 7), title='Preliminary heat drop distribution'):
        logger.info('%s  plot_heat_drop_distribution' % (self.str()))
        x_arr = list(range(1, self.stage_number + 1))
        y_arr = [item.H0 for item in self._stages]
        plt.figure(figsize=figsize)
        plt.plot(x_arr, y_arr, 'o', color='red', markersize=12)
        plt.plot(x_arr, y_arr, color='red')
        plt.grid()
        plt.xticks(x_arr, x_arr)
        plt.xlim(0.7, self.stage_number + 0.3)
        plt.ylim(0, max(y_arr) + 2e4)
        plt.ylabel(r'$H_i$', fontsize=20)
        plt.xlabel(r'$Stage\ number$', fontsize=20)
        plt.title(title, fontsize=20)
        plt.show()

    def _compute_n_max(self):
        logger.info('%s  _compute_n_max' % (self.str()))
        self.n_max = np.sqrt(self.sigma_l / (self.k_n * self.last.A2))
        self.n = self.n_rel * self.n_max
        for item in self._stages:
            item.n = self.n

    def _compute_outlet_static_parameters(self):
        logger.info('%s  _compute_outlet_static_parameters' % (self.str()))
        self.work_fluid.__init__()
        self.work_fluid.alpha = self.alpha_air
        self.rho_t_stag = self.p_t_stag / (self.work_fluid.R * self.T_t_stag)
        self.work_fluid.T1 = self.T_t_stag
        self.dT_t_rel = 1
        self._iter_number_static_par = 0
        while self.dT_t_rel >= 0.001:
            self._iter_number_static_par += 1
            logger.debug('%s _compute_outlet_static_parameters _iter_number = %s' %
                         (self.str(), self._iter_number_static_par))
            self.c_p_gas_t = self.work_fluid.c_p_av_int
            self.k_gas_t = self.work_fluid.k_av_int
            logger.debug('%s _compute_outlet_static_parameters k_gas_t = %s' % (self.str(), self.k_gas_t))
            self.a_cr_t = GasDynamicFunctions.a_cr(self.T_t_stag, self.work_fluid.k_av_int, self.work_fluid.R)

            def eps(c):
                return GasDynamicFunctions.eps_lam(c / self.a_cr_t, self.work_fluid.k)

            def func_to_solve(x):
                return [x[0] * eps(x[0]) - self.G_turbine / (self.rho_t_stag * self.last.A2)]

            x = fsolve(func_to_solve, np.array([200]))
            self.c_t = x[0]
            logger.debug('%s _compute_outlet_static_parameters c_t = %s' % (self.str(), self.c_t))
            self.lam_t = self.c_t / self.a_cr_t
            self.p_t = self.p_t_stag * GasDynamicFunctions.pi_lam(self.lam_t, self.work_fluid.k_av_int)
            self.T_t = self.T_t_stag * GasDynamicFunctions.tau_lam(self.lam_t, self.work_fluid.k_av_int)
            logger.debug('%s _compute_outlet_static_parameters T_t = %s' % (self.str(), self.T_t))
            self.dT_t_rel = abs(self.T_t - self.work_fluid.T2) / self.work_fluid.T2
            logger.debug('%s _compute_outlet_static_parameters dT_t_rel = %s' % (self.str(), self.dT_t_rel))
            self.work_fluid.T2 = self.T_t
        self.work_fluid.T1 = self.T_g_stag
        self.work_fluid.T2 = self.T_t
        self.k_gas = self.work_fluid.k_av_int
        self.c_p_gas = self.work_fluid.c_p_av_int
        self.H_t = self.c_p_gas * self.T_g_stag * (1 - (self.p_g_stag / self.p_t) ** ((1 - self.k_gas) / self.k_gas))
        self.eta_l = func.eta_turb_l(self.eta_t_stag, self.H_t_stag, self.H_t, self.c_t)
        self.alpha = (self.stage_number - 1) / (2 * self.stage_number) * (1 - self.eta_l) * \
                     ((self.p_g_stag / self.p_t) ** ((self.k_gas - 1) / self.k_gas) - 1)

    def compute_output(self, compute_heat_drop_auto=True):
        logger.info('%s  compute_output' % (self.str()))
        self._compute_geometry()
        self._compute_n_max()
        self._compute_outlet_static_parameters()
        if compute_heat_drop_auto is True:
            self._compute_heat_drop_distribution()
        else:
            pass


def specify_h01(turbine_geometry: TurbineGeomAndHeatDropDistribution):
    logger.info('make_h01_more_accurate')
    dh01_rel = 1
    H01 = turbine_geometry.H01
    iter_number = 0
    while dh01_rel >= 0.01:
        iter_number += 1
        logger.debug('make_h01_more_accurate iter_number = %s' % iter_number)
        turbine_geometry.H01 = H01
        logger.debug('make_h01_more_accurate H01 = %s' % H01)
        turbine_geometry.compute_output()
        dh01_rel = abs(turbine_geometry.first.H0 - turbine_geometry.H01) / turbine_geometry.H01
        logger.debug('make_h01_more_accurate dh01_rel = %s' % dh01_rel)
        H01 = turbine_geometry.first.H0


if __name__ == '__main__':
    deg = np.pi / 180
    # stage_geom = _StageGeomAndHeatDrop()
    # stage_geom.D1 = 0.40
    # stage_geom.delta_a_b_rk_ratio = 0.25
    # stage_geom.delta_a_b_sa_ratio = 0.25
    # stage_geom.delta_r_rk = 0.002
    # stage_geom.l1_b_sa_ratio = 2.2
    # stage_geom.l2_b_rk_ratio = 2.5
    # stage_geom.l1 = 0.07
    # stage_geom.gamma_av = 5 * deg
    # stage_geom.gamma_out = 8 * deg
    # stage_geom.x0 = 0.02
    # stage_geom.gamma_in = np.arctan(np.tan(stage_geom.gamma_out) - 2 * np.tan(stage_geom.gamma_av))
    # stage_geom.plot()
    # plt.grid()
    # plt.show()
    # stage_geom.n_rel = 20e3
    turbine_geom = TurbineGeomAndHeatDropDistribution(2, 0.91, 300e3, 150, 0.9, KeroseneCombustionProducts(),
                                                      1400, 1.3e6, 2.4, 9, 1 / 3.5, 180e3, 0.8, 0.96, 15 * deg,
                                                      6.8, 200e6, 850, 120e3, gamma_av=0 * deg,
                                                      gamma_sum=15 * deg)
    turbine_geom[0].delta_a_b_sa_ratio = 0.23
    turbine_geom[0].delta_a_b_rk_ratio = 0.25
    turbine_geom[0].l1_b_sa_ratio = 2
    turbine_geom[0].l2_b_rk_ratio = 2.5
    turbine_geom[0].delta_r_rk_l2_ratio = 0.01
    turbine_geom[0].mu = 1
    turbine_geom[1].delta_a_b_sa_ratio = 0.25
    turbine_geom[1].delta_a_b_rk_ratio = 0.28
    turbine_geom[1].delta_r_rk_l2_ratio = 0.01
    turbine_geom[1].l1_b_sa_ratio = 2.5
    turbine_geom[1].l2_b_rk_ratio = 2.8
    turbine_geom[1].mu = 1
    specify_h01(turbine_geom)
    turbine_geom.plot_geometry()
    turbine_geom.plot_heat_drop_distribution(figsize=(10, 7))
    print(turbine_geom.c_t)
    print(turbine_geom.p_t)
    print(turbine_geom.T_t)
    print(turbine_geom.lam_t)
    print(turbine_geom.n)
    print(turbine_geom.gamma_av / deg)
    print(turbine_geom.gamma_out / deg)
    print(turbine_geom.gamma_in / deg)
    print(turbine_geom.D1)
    print(turbine_geom.l1)
    print(turbine_geom[1].D1)
    print(turbine_geom[1].l1)
    print(turbine_geom.alpha)





