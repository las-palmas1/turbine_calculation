from average_streamline_calculation.turbine_stage_geometry import InvalidStageSizeValue, StageGeomAndHeatDrop, \
    TurbineGeomAndHeatDropDistribution
import functions as func
import logging
from gas import KeroseneCombustionProducts
import numpy as np
import matplotlib.pyplot as plt
import os

log_filename = os.path.join(os.path.dirname(__file__), 'average_streamline_calculation.log')
logger = func.create_logger(__name__, logging.DEBUG, filename=log_filename, filemode='a',
                            add_console_handler=False)


class StageGasDynamics:
    def __init__(self, T0_stag, p0_stag, G_stage_in, G_turbine, alpha_air, work_fluid: KeroseneCombustionProducts, rho,
                 phi, psi, l1, l2, D1, D2, delta_r_rk, n, epsilon, g_lk, g_ld, g_lb, **kwargs):
        """

        :param T0_stag:
        :param p0_stag:
        :param G_stage_in: расход газа на входе в ступень
        :param G_turbine: расход газа через СА первой ступени
        :param alpha_air:
        :param work_fluid:
        :param rho:
        :param phi:
        :param psi:
        :param l1:
        :param l2:
        :param D1:
        :param D2:
        :param delta_r_rk:
        :param n:
        :param epsilon:
        :param g_lk: относительный расход утечек в концевых лабиринтах
        :param g_ld: относительный расход перетечек в лабиринтных уплотнениях сопловых диафрагм
        :param g_lb: относительный расход перетечек поверх бондажа рабочих лопаток
        :param kwargs: H0, p2, L_t, eta_t0
        """
        self._T0_stag = T0_stag
        self._p0_stag = p0_stag
        self._G_stage_in = G_stage_in
        self._G_turbine = G_turbine
        self._work_fluid = work_fluid
        self._alpha_air = alpha_air
        self._phi = phi
        self._psi = psi
        self._l1 = l1
        self._l2 = l2
        self._D1 = D1
        self._D2 = D2
        self._rho = rho
        self._n = n
        self._delta_r_rk = delta_r_rk
        self._epsilon = epsilon
        self._g_lk = g_lk
        self._g_ld = g_ld
        self._g_lb = g_lb
        self._kwargs = kwargs
        if 'H0' in kwargs:
            self._H0 = kwargs['H0']
            self._p2 = None
            self._L_t = None
            self._eta_t0 = None
        elif 'p2' in kwargs:
            self._H0 = None
            self._p2 = kwargs['p2']
            self._L_t = None
            self._eta_t0 = None
        elif ('L_t' in kwargs) and ('eta_t0' in kwargs):
            self._L_t = kwargs['L_t']
            self._eta_t0 = kwargs['eta_t0']
            self._H0 = None
            self._p2 = None
        else:
            assert False, 'H0 or p2 or (L_t and eta_t0) must be set'
        self._stage_calculation()

    @property
    def G_turbine(self):
        return self._G_turbine

    @G_turbine.setter
    def G_turbine(self, value):
        self._G_turbine = value
        self._stage_calculation()

    @property
    def g_lb(self):
        return self._g_lb

    @g_lb.setter
    def g_lb(self, value):
        self._g_lb = value
        self._stage_calculation()

    @property
    def g_ld(self):
        return self._g_ld

    @g_ld.setter
    def g_ld(self, value):
        self._g_ld = value
        self._stage_calculation()

    @property
    def g_lk(self):
        return self._g_lk

    @g_lk.setter
    def g_lk(self, value):
        self._g_lk = value
        self._stage_calculation()

    @property
    def eta_t0(self):
        if ('L_t' in self._kwargs) and ('eta_t0' in self._kwargs):
            assert self._eta_t0 is not None, 'eta_t0 must not be None'
        return self._eta_t0

    @eta_t0.setter
    def eta_t0(self, value):
        self._eta_t0 = value
        if ('L_t' in self._kwargs) and ('eta_t0' in self._kwargs):
            self._stage_calculation()

    @property
    def L_t(self):
        if ('L_t' in self._kwargs) and ('eta_t0' in self._kwargs):
            assert self._L_t is not None, 'L_t must not be None'
        return self._L_t

    @L_t.setter
    def L_t(self, value):
        self._L_t = value
        if ('L_t' in self._kwargs) and ('eta_t0' in self._kwargs):
            self._stage_calculation()

    @property
    def H0(self):
        if 'H0' in self._kwargs:
            assert self._H0 is not None, 'H0 must not be None'
        return self._H0

    @H0.setter
    def H0(self, value):
        self._H0 = value
        if 'H0' in self._kwargs:
            self._stage_calculation()

    @property
    def p2(self):
        if 'p2' in self._kwargs:
            assert self._p2 is not None, 'p2 must not be None'
        return self._p2

    @p2.setter
    def p2(self, value):
        self._p2 = value
        if 'p2' in self._kwargs:
            self._stage_calculation()

    @property
    def work_fluid(self) -> KeroseneCombustionProducts:
        assert self._work_fluid is not None, 'work_fluid can not be None'
        return self._work_fluid

    @property
    def epsilon(self):
        return self._epsilon

    @epsilon.setter
    def epsilon(self, value):
        self._epsilon = value
        self._stage_calculation()

    @property
    def delta_r_rk(self):
        return self._delta_r_rk

    @delta_r_rk.setter
    def delta_r_rk(self, value):
        self._delta_r_rk = value
        self._stage_calculation()

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, value):
        self._n = value
        self._stage_calculation()

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, value):
        self._rho = value
        self._stage_calculation()

    @property
    def D2(self):
        return self._D2

    @D2.setter
    def D2(self, value):
        self._D2 = value
        self._stage_calculation()

    @property
    def D1(self):
        return self._D1

    @D1.setter
    def D1(self, value):
        self._D1 = value
        self._stage_calculation()

    @property
    def l2(self):
        return self._l2

    @l2.setter
    def l2(self, value):
        self._l2 = value
        self._stage_calculation()

    @property
    def l1(self):
        return self._l1

    @l1.setter
    def l1(self, value):
        self._l1 = value
        self._stage_calculation()

    @property
    def psi(self):
        return self._psi

    @psi.setter
    def psi(self, value):
        self._psi = value
        self._stage_calculation()

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, value):
        self._phi = value
        self._stage_calculation()

    @property
    def alpha_air(self):
        return self._alpha_air

    @alpha_air.setter
    def alpha_air(self, value):
        self._alpha_air = value
        self._stage_calculation()

    @property
    def G_stage_in(self):
        return self._G_stage_in

    @G_stage_in.setter
    def G_stage_in(self, value):
        self._G_stage_in = value
        self._stage_calculation()

    @property
    def T0_stag(self):
        return self._T0_stag

    @T0_stag.setter
    def T0_stag(self, value):
        self._T0_stag = value
        self._stage_calculation()

    @property
    def p0_stag(self):
        return self._p0_stag

    @p0_stag.setter
    def p0_stag(self, value):
        self._p0_stag = value
        self._stage_calculation()

    def str(self):
        str_arr = str(self).split()
        return str_arr[0][1:] + ' ' + str_arr[1]

    def _stage_calculation(self):
        logger.info('%s _stage_calculation' % self.str())
        if 'H0' in self._kwargs:
            self._specified_heat_drop_calculation()
        elif 'p2' in self._kwargs:
            self._specified_outlet_pressure_calculation()
        elif 'L_t' in self._kwargs and 'eta_t0' in self._kwargs:
            self._specified_work_calculation()

    def _specified_heat_drop_calculation(self):
        """Вычисляет параметры ступени по заданному теплоперепаду с уточнением k"""
        logger.info('%s _specified_heat_drop_calculation' % self.str())
        self.work_fluid.__init__()
        self.work_fluid.T1 = self.T0_stag
        self.work_fluid.alpha = self.alpha_air
        self.dk_rel = 1
        self._iter_number_k_gas = 0
        while self.dk_rel >= 0.001:
            self._iter_number_k_gas += 1
            logger.debug('%s _specified_heat_drop_calculation _iter_number_k_gas = %s' %
                         (self.str(), self._iter_number_k_gas))
            self._compute_stage_parameters()

    def _specified_outlet_pressure_calculation(self):
        """Вычисляет параметры ступени по заданному выходному давлению"""
        logger.info('%s _specified_outlet_pressure_calculation' % self.str())
        self.work_fluid.__init__()
        self.work_fluid.T1 = self.T0_stag
        self.work_fluid.alpha = self.alpha_air
        self.dk_rel = 1
        self._iter_number_k_gas = 0
        while self.dk_rel >= 0.001:
            self._iter_number_k_gas += 1
            logger.debug('%s _specified_outlet_pressure_calculation _iter_number_k_gas = %s' %
                         (self.str(), self._iter_number_k_gas))
            self.H0 = self.work_fluid.c_p_av_int * \
                       self.T0_stag * (1 - (self.p0_stag / self.p2) **
                                       ((1 - self.work_fluid.k_av_int) / self.work_fluid.k_av_int))
            logger.debug('%s _specified_outlet_pressure_calculation H0 = %s' % (self.str(), self.H0))
            self._compute_stage_parameters()

    def _specified_work_calculation(self):
        """Вычисляет параметры ступени по заданной работе"""
        logger.info('%s _specified_work_calculation' % self.str())
        self.d_eta_t_rel = 1
        eta_t = self._eta_t0
        self._iter_number_eta_t = 0
        while self.d_eta_t_rel >= 0.001:
            self._iter_number_eta_t += 1
            logger.debug('%s _specified_work_calculation _iter_number_eta_t = %s' %
                         (self.str(), self._iter_number_eta_t))
            self.H0 = self.L_t / eta_t
            self._specified_heat_drop_calculation()
            self.d_eta_t_rel = abs(self.eta_t - eta_t) / eta_t
            eta_t = self.eta_t

    def _compute_stage_parameters(self):
        """Вычисляет параметры ступени по известному теплоперепаду без уточнения k"""
        logger.info('%s _compute_stage_parameters' % self.str())
        self.u1 = np.pi * self.D1 * self.n / 60
        self.H_s = self.H0 * (1 - self.rho)
        self.c1 = self.phi * np.sqrt(2 * self.H_s)
        self.c_p_gas = self.work_fluid.c_p_av_int
        self.T1 = self.T0_stag - self.H_s * self.phi ** 2 / self.c_p_gas
        self.T1_ad = self.T0_stag - self.H_s / self.c_p_gas
        self.k_gas = self.work_fluid.k_av_int
        logger.debug('%s _compute_stage_parameters k = %s' % (self.str(), self.k_gas))
        self.p1 = self.p0_stag * (self.T1_ad / self.T0_stag) ** (self.k_gas / (self.k_gas - 1))
        self.A1_a = np.pi * self.D1 * self.l1
        self.rho1 = self.p1 / (self.work_fluid.R * self.T1)
        self.G_sa = self.G_stage_in - self.g_ld * self.G_turbine
        self.c1_a = self.G_sa / (self.rho1 * self.A1_a)
        if self.c1_a > self.c1:
            raise InvalidStageSizeValue('c1_a must be less than c1')
        self.alpha1 = np.arcsin(self.c1_a / self.c1)
        self.c1_u = self.c1 * np.cos(self.alpha1)
        self.w1 = np.sqrt(self.c1**2 + self.u1**2 - 2 * self.c1 * self.u1 * np.cos(self.alpha1))
        if self.c1 * np.cos(self.alpha1) - self.u1 >= 0:
            self.beta1 = np.arctan(self.c1_a / (self.c1 * np.cos(self.alpha1) - self.u1))
        else:
            self.beta1 = np.pi + np.arctan(self.c1_a / (self.c1 * np.cos(self.alpha1) - self.u1))
        self.H_l = self.rho * self.H0 * self.T1 / self.T1_ad
        self.u2 = np.pi * self.D2 * self.n / 60
        self.w2 = self.psi * np.sqrt(self.w1**2 + 2 * self.H_l + self.u2**2 - self.u1**2)
        self.T2 = self.T1 - (self.w2**2 - self.w1**2 - self.u2**2 + self.u1**2) / (2 * self.c_p_gas)
        self.T2_ad = self.T1 - self.H_l / self.c_p_gas
        if ('H0' in self._kwargs) or ('L_t' in self._kwargs and 'eta_t0' in self._kwargs):
            self.p2 = self.p1 * (self.T2_ad / self.T1) ** (self.k_gas / (self.k_gas - 1))
        elif 'p2' in self._kwargs:
            self.p2_check = self.p1 * (self.T2_ad / self.T1) ** (self.k_gas / (self.k_gas - 1))
        self.rho2 = self.p2 / (self.work_fluid.R * self.T2)
        self.A2_a = np.pi * self.D2 * self.l2
        self.G_rk = self.G_stage_in - self.G_turbine * (self.g_lb + self.g_lk)
        self.c2_a = self.G_rk / (self.A2_a * self.rho2)
        if self.c2_a / self.w2 >= 1:
            raise InvalidStageSizeValue('c2_a must be less than w2')
        self.beta2 = np.arcsin(self.c2_a / self.w2)
        self.c2_u = self.w2 * np.cos(self.beta2) - self.u2
        if self.c2_u >= 0:
            self.alpha2 = np.arctan(self.c2_a / self.c2_u)
        else:
            self.alpha2 = np.pi + np.arctan(self.c2_a / self.c2_u)
        self.c2 = np.sqrt(self.c2_a**2 + self.c2_u**2)
        self.L_u = self.c1_u * self.u1 + self.c2_u * self.u2
        self.eta_u = self.L_u / self.H0
        self.h_s = (1 / self.phi ** 2 - 1) * self.c1 ** 2 / 2
        self.h_s_touch = self.h_s * self.T2_ad / self.T1
        self.zeta_s = self.h_s / self.H0
        self.zeta_s_touch = self.h_s_touch / self.H0
        self.h_l = (1 / self.psi ** 2 - 1) * self.w2 ** 2 / 2
        self.zeta_l = self.h_l / self.H0
        self.h_v = self.c2 ** 2 / 2
        self.zeta_v = self.h_v / self.H0
        self.eta_u_check = 1 - self.zeta_s_touch - self.zeta_l - self.zeta_v
        self.D_av = 0.5 * (self.D1 + self.D2)
        self.h_z = 1.37 * (1 + 1.6 * self.rho) * (1 + self.l2 / self.D_av) * self.delta_r_rk / self.l2 * self.L_u
        self.zeta_z = self.h_z / self.H0
        self.L_uz = self.L_u - self.h_z
        self.eta_t_touch = self.eta_u - self.zeta_z
        self.eta_l_touch = self.eta_t_touch + self.zeta_v
        self.l_av = 0.5 * (self.l1 + self.l2)
        self.u_av = 0.5 * (self.u1 + self.u2)
        self.N_tv = (1.07 * self.D_av**2 + 61 * (1 - self.epsilon) * self.D_av * self.l_av) * \
                    (self.u_av / 100)**3 * self.rho
        self.h_tv = self.N_tv / self.G_stage_in
        self.zeta_tv = self.h_tv / self.H0
        self.eta_t = self.eta_t_touch - self.zeta_tv
        self.eta_l = self.eta_l_touch - self.zeta_tv
        if ('H0' in self._kwargs) or ('p2' in self._kwargs):
            self.L_t = self.H0 * self.eta_t
            logger.debug('%s _compute_stage_parameters L_t = %s' % (self.str(), self.L_t))
            # удельная работа ступени, отнесенная к расходу через СА первой ступени с учетом потерь из-за утечек
            self.L_t_rel = self.L_t * self.G_stage_in / self.G_turbine - self.L_t * (self.g_ld + self.g_lk + self.g_lb)
        self.T_st = self.T2 + self.h_z / self.c_p_gas + self.h_tv / self.c_p_gas
        self.T_st_stag = self.T_st + self.h_v / self.c_p_gas
        self.work_fluid.T2 = self.T_st
        try:
            self.dk_rel = abs(self.k_gas - self.work_fluid.k_av_int) / self.k_gas
        except TypeError:
            self.dk_rel = max(abs(self.k_gas - self.work_fluid.k_av_int) / self.k_gas)
        logger.debug('%s _compute_stage_parameters dk_rel = %s' % (self.str(), self.dk_rel))
        self.p2_stag = self.p2 * (self.T_st_stag / self.T_st) ** (self.k_gas / (self.k_gas - 1))
        self.H0_stag = self.c_p_gas * self.T0_stag * (1 - (self.p2_stag / self.p2) ** ((self.k_gas - 1) / self.k_gas))
        self.eta_t_stag = self.L_t / self.H0_stag
        self.G_stage_out = self.G_stage_in - self.G_turbine * self.g_lk

    def plot_velocity_triangle(self, title='', figsize=(8, 8)):
        x_in = np.array([0, -self.c1_u, -self.c1_u + self.u1, 0])
        y_in = np.array([self.c1_a, 0, 0, self.c1_a])
        x_out = np.array([0, self.c2_u, self.c2_u + self.u2, 0])
        y_out = np.array([self.c1_a, self.c1_a - self.c2_a, self.c1_a - self.c2_a, self.c1_a])
        plt.figure(figsize=figsize)
        plt.plot(x_in, y_in, linewidth=2, color='red', label='inlet')
        plt.plot(x_out, y_out, linewidth=2, color='blue', label='outlet')
        plt.xlim(-self.c1_u, self.c2_u + self.u2)
        plt.ylim(-max(self.c1_a, self.c1_u), max(self.c1_a, self.c2_u + self.u2))
        plt.grid()
        plt.title(title, fontsize=20)
        plt.legend()
        plt.show()


def get_first_stage_gas_dynamics(stage_geom: StageGeomAndHeatDrop, T0_stag, p0_stag, G_turbine,
                                 alpha_air) -> StageGasDynamics:
    logger.info('get_first_stage_gas_dynamics')
    result = StageGasDynamics(T0_stag, p0_stag, G_turbine, G_turbine, alpha_air, KeroseneCombustionProducts(),
                              stage_geom.rho, stage_geom.phi, stage_geom.psi, stage_geom.l1,
                              stage_geom.l2, stage_geom.D1, stage_geom.D2, stage_geom.delta_r_rk,
                              stage_geom.n, stage_geom.epsilon, stage_geom.g_lk, stage_geom.g_ld,
                              stage_geom.g_lb, H0=stage_geom.H0)
    return result


def get_intermediate_stage(stage_geom: StageGeomAndHeatDrop, prev_stage: StageGasDynamics, precise_heat_drop=True) -> \
        StageGasDynamics:
    logger.info('get_intermediate_stage')
    if precise_heat_drop:
        H0 = stage_geom.H0 * (1 + (1 - stage_geom.mu) ** 2 * prev_stage.c2 ** 2 /
                              (2 * prev_stage.c_p_gas * prev_stage.T_st)) + 0.5 * (stage_geom.mu * prev_stage.c2) ** 2
    else:
        H0 = stage_geom.H0
    p0_stag = prev_stage.p2 * (1 + (stage_geom.mu * prev_stage.c2)**2 / (2 * prev_stage.c_p_gas * prev_stage.T_st)) ** \
                              (prev_stage.k_gas / (prev_stage.k_gas - 1))
    result = StageGasDynamics(prev_stage.T_st_stag, p0_stag, prev_stage.G_stage_out, prev_stage.G_turbine,
                              prev_stage.alpha_air, KeroseneCombustionProducts(), stage_geom.rho, stage_geom.phi,
                              stage_geom.psi, stage_geom.l1, stage_geom.l2, stage_geom.D1, stage_geom.D2,
                              stage_geom.delta_r_rk, stage_geom.n, stage_geom.epsilon, stage_geom.g_lk,
                              stage_geom.g_ld, stage_geom.g_lb, H0=H0)
    return result


def get_last_pressure_stage(turbine_geom: TurbineGeomAndHeatDropDistribution,
                            stage_geom: StageGeomAndHeatDrop, prev_stage: StageGasDynamics) -> StageGasDynamics:
    logger.info('get_last_pressure_stage')
    p0_stag = prev_stage.p2 * (1 + (stage_geom.mu * prev_stage.c2)**2 / (2 * prev_stage.c_p_gas * prev_stage.T_st)) ** \
                              (prev_stage.k_gas / (prev_stage.k_gas - 1))
    result = StageGasDynamics(prev_stage.T_st_stag, p0_stag, prev_stage.G_stage_out, prev_stage.G_turbine,
                              prev_stage.alpha_air, KeroseneCombustionProducts(), stage_geom.rho, stage_geom.phi,
                              stage_geom.psi, stage_geom.l1, stage_geom.l2, stage_geom.D1, stage_geom.D2,
                              stage_geom.delta_r_rk, stage_geom.n, stage_geom.epsilon, stage_geom.g_lk,
                              stage_geom.g_ld, stage_geom.g_lb, p2=turbine_geom.p_t)
    return result


def get_last_work_stage(stage_geom: StageGeomAndHeatDrop, prev_stage: StageGasDynamics,
                        L_stage_rel, eta_t0) -> StageGasDynamics:
    logger.info('get_last_work_stage')
    L_stage = L_stage_rel / (prev_stage.G_stage_out / prev_stage.G_turbine -
                             (stage_geom.g_lb + stage_geom.g_ld + stage_geom.g_lk))
    p0_stag = prev_stage.p2 * (1 + (stage_geom.mu * prev_stage.c2) ** 2 / (2 * prev_stage.c_p_gas *
                                                                           prev_stage.T_st)) ** \
                              (prev_stage.k_gas / (prev_stage.k_gas - 1))
    result = StageGasDynamics(prev_stage.T_st_stag, p0_stag, prev_stage.G_stage_out, prev_stage.G_turbine,
                              prev_stage.alpha_air, KeroseneCombustionProducts(), stage_geom.rho, stage_geom.phi,
                              stage_geom.psi, stage_geom.l1, stage_geom.l2, stage_geom.D1, stage_geom.D2,
                              stage_geom.delta_r_rk, stage_geom.n, stage_geom.epsilon, stage_geom.g_lk,
                              stage_geom.g_ld, stage_geom.g_lb, L_t=L_stage, eta_t0=eta_t0)
    return result


if __name__ == '__main__':
    deg = np.pi / 180
    stage_gd = StageGasDynamics(1400, 1100e3, 30, 30, 2.5, KeroseneCombustionProducts(), 0.8, 0.96, 0.96, 0.08, 0.10,
                                0.34, 0.37, 0.001, 15e3, 1, 0, 0, 0, L_t=250e3, eta_t0=0.85)
    print(stage_gd.alpha1 / deg)
    print(stage_gd.alpha2 / deg)
    print(stage_gd.beta1 / deg)
    print(stage_gd.beta2 / deg)
    print(stage_gd.H_l)
    print(stage_gd.H_s)
    print(stage_gd.c1_a)
    print(stage_gd.c2_a)
    print(stage_gd.T2)
    print(stage_gd.T_st)
    print(stage_gd.T_st_stag)
    print(stage_gd.eta_t)
    stage_gd.plot_velocity_triangle()
    stage_gd.L_t = 190e3
    print(stage_gd.eta_t)