import logging.config

import pandas

import functions as func
import gas_dynamics as gd
from cycle_calculation.gas_turbine_units_input import *
from cycle_calculation.gas_turbine_units_input import OutletInput
from cycle_calculation.gas_turbine_units_output import *
import os

log_file_name = os.path.join(os.path.dirname(__file__), 'cycle_calculation.log')
logger = func.create_logger(loggerlevel=logging.INFO, filename=log_file_name, filemode='w', name=__name__)


def compute_comp_output(work_fluid: Air, pi_comp_stag: np.array, eta_comp_stag_p, T_in_stag):
    work_fluid.__init__()
    dk_rel = 1
    k_air_ret = np.array([work_fluid.k])
    work_fluid.T1 = T_in_stag
    eta_comp_stag = None
    k = 0
    k_air_old = None
    while dk_rel >= 0.01:
        k += 1
        eta_comp_stag = func.eta_comp_stag(pi_comp_stag, k_air_ret, eta_comp_stag_p)
        work_fluid.T2 = T_in_stag * (1 + (pi_comp_stag ** ((k_air_ret - 1) / k_air_ret) - 1) / eta_comp_stag)
        k_air_ret_new = work_fluid.k_av_int
        dk_rel = max(abs(k_air_ret - k_air_ret_new) / k_air_ret)
        k_air_old = k_air_ret
        k_air_ret = k_air_ret_new
    L_comp = work_fluid.c_p_av_int * (work_fluid.T2 - work_fluid.T1)
    return {'k_air': k_air_ret, 'c_p_air': work_fluid.c_p_av_int, 'eta_comp_stag': eta_comp_stag,
            'T_out_stag': work_fluid.T2, 'L_comp': L_comp, 'work_fluid': work_fluid, 'k_air_old': k_air_old,
            'dk_rel': dk_rel}


class Compressor:
    def __init__(self):
        self._input = CompressorInput()
        self._output = CompressorOutput()

    def _compute_comp_output(self, work_fluid: Air, pi_comp_stag: np.array, eta_comp_stag_p, T_in_stag):
        self.work_fluid = work_fluid
        self.work_fluid.__init__()
        self.pi_comp_stag = pi_comp_stag
        self.eta_comp_stag_p = eta_comp_stag_p
        self.T_in_stag = T_in_stag
        self.dk_rel = 1
        self.k_air_ret = np.array([work_fluid.k])
        self.work_fluid.T1 = T_in_stag
        self.eta_comp_stag = None
        self.k = 0
        self.k_air_old = None
        while self.dk_rel >= 0.01:
            self.k += 1
            self.eta_comp_stag = func.eta_comp_stag(self.pi_comp_stag, self.k_air_ret, self.eta_comp_stag_p)
            work_fluid.T2 = self.T_in_stag * (1 + (self.pi_comp_stag ** ((self.k_air_ret - 1) / self.k_air_ret) - 1) /
                                              self.eta_comp_stag)
            self.k_air_ret_new = self.work_fluid.k_av_int
            self.dk_rel = max(abs(self.k_air_ret - self.k_air_ret_new) / self.k_air_ret)
            self.k_air_old = self.k_air_ret
            self.k_air_ret = self.k_air_ret_new
        self.L_comp = self.work_fluid.c_p_av_int * (self.work_fluid.T2 - self.work_fluid.T1)
        return {'k_air': self.k_air_ret, 'c_p_air': self.work_fluid.c_p_av_int, 'eta_comp_stag': self.eta_comp_stag,
                'T_out_stag': self.work_fluid.T2, 'L_comp': self.L_comp, 'work_fluid': self.work_fluid,
                'k_air_old': self.k_air_old,
                'dk_rel': self.dk_rel}

    @property
    def input(self) -> CompressorInput:
        return self._input

    @input.setter
    def input(self, value: CompressorInput):
        self._input = value

    @property
    def output(self) -> CompressorOutput:
        return self._output

    def compute_output(self):
        output = self._compute_comp_output(self._input.work_fluid, self._input.pi_comp_stag,
                                           self._input.eta_comp_stag_p, self._input.T_in_stag)
        self._output.k_air = self.k_air_ret
        self._output.p_out_stag = self._input.p_in_stag * self._input.pi_comp_stag
        self._output.c_p_air = self.work_fluid.c_p_av_int
        self._output.eta_comp_stag = self.eta_comp_stag
        self._output.T_out_stag = self.work_fluid.T2
        self._output.L_comp = self.L_comp
        self._output.work_fluid = self.work_fluid
        output['p_out_stag'] = self._output.p_out_stag
        del(output['work_fluid'])
        self._output.output_frame = pandas.DataFrame.from_dict(output)


def compute_power_turb_output(work_fluid: KeroseneCombustionProducts, T_out_stag_init, T_in_stag, p_in_stag, p_out_stag,
                              eta_turb_stag_p, alpha):
    work_fluid.__init__()
    pi_turb_stag = p_in_stag / p_out_stag
    dT_rel = 1
    work_fluid.T1 = T_in_stag
    work_fluid.T2 = np.array([T_out_stag_init])
    work_fluid.alpha = alpha
    eta_turb_stag = None
    while dT_rel >= 0.01:
        eta_turb_stag = func.eta_turb_stag(pi_turb_stag, work_fluid.k_av_int, eta_turb_stag_p)
        k_air_old = work_fluid.k_av_int
        work_fluid.T2 = T_in_stag * (1 - (1 - pi_turb_stag **
                                          ((1 - work_fluid.k_av_int) / work_fluid.k_av_int)) * eta_turb_stag)
        k_air_new = work_fluid.k_av_int
        dT_rel = max(abs(work_fluid.T2 - work_fluid.T2) / work_fluid.T2)
    L_turb = work_fluid.c_p_av_int * (T_in_stag - work_fluid.T2)
    H_turb_stag = work_fluid.c_p_av_int * \
                  (1 - pi_turb_stag**((1 - work_fluid.k_av_int) / work_fluid.k_av_int)) * T_in_stag
    return {'k': work_fluid.k_av_int, 'c_p': work_fluid.c_p_av_int, 'pi_turb_stag': pi_turb_stag,
            'T_out_stag': work_fluid.T2, 'L_turb': L_turb, 'H_turb_stag': H_turb_stag,
            'eta_turb_stag': eta_turb_stag, 'work_fluid': work_fluid}


class PowerTurbine:
    def __init__(self):
        self._input = PowerTurbineInput()
        self._output = PowerTurbineOutput()

    def _compute_power_turb_output(self, work_fluid: KeroseneCombustionProducts, T_out_stag_init, T_in_stag, p_in_stag,
                                   p_out_stag, eta_turb_stag_p, alpha):
        self.work_fluid = work_fluid
        self.work_fluid.__init__()
        self.T_out_stag_init = T_out_stag_init
        self.T_in_stag = T_in_stag
        self.p_in_stag = p_in_stag
        self.p_out_stag = p_out_stag
        self.eta_turb_stag_p = eta_turb_stag_p
        self.alpha = alpha
        self.pi_turb_stag = p_in_stag / p_out_stag
        self.dT_rel = 1
        self.work_fluid.T1 = T_in_stag
        self.work_fluid.T2 = np.array([T_out_stag_init])
        self.work_fluid.alpha = alpha
        self.eta_turb_stag = None
        while self.dT_rel >= 0.01:
            self.eta_turb_stag = func.eta_turb_stag(self.pi_turb_stag, self.work_fluid.k_av_int, self.eta_turb_stag_p)
            self.k_gas_old = self.work_fluid.k_av_int
            self.c_p_gas_old = self.work_fluid.c_p_av_int
            self.work_fluid.T2 = self.T_in_stag * (1 - (1 - self.pi_turb_stag **
                                                   ((1 - self.work_fluid.k_av_int) / self.work_fluid.k_av_int)) *
                                                   self.eta_turb_stag)
            self.k_gas_new = self.work_fluid.k_av_int
            self.dT_rel = max(abs(self.work_fluid.T2 - self.work_fluid.T2) / self.work_fluid.T2)
            self.dk_rel = max(abs(self.k_gas_new - self.k_gas_old) / self.k_gas_old)
        self.L_turb = self.work_fluid.c_p_av_int * (self.T_in_stag - self.work_fluid.T2)
        self.H_turb_stag = self.work_fluid.c_p_av_int * \
                      (1 - self.pi_turb_stag ** ((1 - self.work_fluid.k_av_int) / self.work_fluid.k_av_int)) * \
                           self.T_in_stag
        return {'k': self.work_fluid.k_av_int, 'c_p': self.work_fluid.c_p_av_int, 'pi_turb_stag': self.pi_turb_stag,
                'T_out_stag': self.work_fluid.T2, 'L_turb': self.L_turb, 'H_turb_stag': self.H_turb_stag,
                'eta_turb_stag': self.eta_turb_stag, 'work_fluid': self.work_fluid}

    @property
    def input(self) -> PowerTurbineInput:
        return self._input

    @input.setter
    def input(self, value: PowerTurbineInput):
        self._input = value

    @property
    def output(self) -> PowerTurbineOutput:
        return self._output

    @staticmethod
    def eta_turb_l(eta_turb_stag, H_turb_stag, H_turb, c_out):
        """Лопаточный КПД турбины через КПД по параметрам торможения"""
        return eta_turb_stag * H_turb_stag / H_turb + c_out ** 2 / (2 * H_turb)

    @staticmethod
    def eta_turb_stag(eta_turb_l, H_turb_stag, H_turb, c_out):
        """КПД по парметрам торможения через лопаточный КПД"""
        return (eta_turb_l * H_turb - c_out ** 2 / 2) / H_turb_stag

    def compute_output(self):
        self._output.pi_turb = self._input.p_in_stag / self._input.p_out_stag
        output = self._compute_power_turb_output(self._input.work_fluid, self._input.T_out_stag_init, self._input.T_in_stag,
                                           self._input.p_in_stag, self.input.p_out_stag, self._input.eta_turb_stag_p,
                                           self._input.alpha)
        self._output.k_gas = output['k']
        self._output.c_p_gas = output['c_p']
        self._output.T_out_stag = output['T_out_stag']
        self._output.eta_turb_stag = output['eta_turb_stag']
        self._output.H_turb_stag = output['H_turb_stag']
        self._output.L_turb = output['L_turb']
        self._output.pi_turb_stag = output['pi_turb_stag']
        self._output.work_fluid = output['work_fluid']
        del(output['work_fluid'])
        self._output.output_frame = pandas.DataFrame.from_dict(output)


def compute_comp_turb_output(work_fluid: KeroseneCombustionProducts, pi_turb_stag_init, T_in_stag, p_in_stag, L_comp,
                             g_gas, eta_turb_stag_p, eta_m, alpha):
    work_fluid.__init__()
    L_turb = L_comp / (eta_m * g_gas)
    logger.debug('comp_turb_output L_turb = %s' % L_turb)
    k_gas = np.array([work_fluid.k])
    dk_rel = 1
    work_fluid.T1 = T_in_stag
    work_fluid.alpha = alpha
    while dk_rel >= 0.01:
        work_fluid.T2 = T_in_stag - L_turb / work_fluid.c_p_av_int
        logger.debug('comp_turb_output T_out_stag = %s' % work_fluid.T2)
        logger.debug('comp_turb_output c_p_gas_av = %s' % work_fluid.c_p_av_int)
        k_gas_new = work_fluid.k_av_int
        dk_rel = max(abs(k_gas_new - k_gas) / k_gas)
        k_gas_old = k_gas
        k_gas = k_gas_new
        logger.debug('comp_turb_output k = %s' % k_gas)
    pi_turb_stag = np.array([pi_turb_stag_init])
    dpi_rel = 1
    p_out_stag = None
    eta_turb_stag = None
    while dpi_rel >= 0.01:
        eta_turb_stag = func.eta_turb_stag(pi_turb_stag, k_gas, eta_turb_stag_p)
        p_out_stag = p_in_stag * (1 -
                                  L_turb / (work_fluid.c_p_av_int * eta_turb_stag * T_in_stag)) ** (k_gas / (k_gas - 1))
        logger.debug('comp_turb_output p_out_stag = %s' % p_out_stag)
        pi_turb_stag_new = p_in_stag / p_out_stag
        dpi_rel = max(abs(pi_turb_stag - pi_turb_stag_new) / pi_turb_stag)
        pi_turb_stag = pi_turb_stag_new
    H_turb_stag = work_fluid.c_p_av_int * T_in_stag * (1 - pi_turb_stag ** ((1 - k_gas) / k_gas))
    return {'L_turb': L_turb, 'k': k_gas, 'c_p': work_fluid.c_p_av_int, 'T_out_stag': work_fluid.T2,
            'p_out_stag': p_out_stag, 'pi_turb_stag': pi_turb_stag, 'H_turb_stag': H_turb_stag,
            'eta_turb_stag': eta_turb_stag, 'work_fluid': work_fluid}


class CompressorTurbine:
    def __init__(self):
        self._input = CompressorTurbineInput()
        self._output = CompressorTurbineOutput()

    def _compute_comp_turb_output(self, work_fluid: KeroseneCombustionProducts, pi_turb_stag_init, T_in_stag, p_in_stag,
                                 L_comp, g_gas, eta_turb_stag_p, eta_m, alpha):
        self.work_fluid = work_fluid
        self.work_fluid.__init__()
        self.pi_turb_stag_init = pi_turb_stag_init
        self.T_in_stag = T_in_stag
        self.p_in_stag = p_in_stag
        self.L_comp = L_comp
        self.g_gas = g_gas
        self.eta_turb_stag_p = eta_turb_stag_p
        self.eta_m = eta_m
        self.alpha = alpha
        self.L_turb = self.L_comp / (self.eta_m * self.g_gas)
        logger.debug('comp_turb_output L_turb = %s' % self.L_turb)
        self.k_gas = np.array([self.work_fluid.k])
        self.dk_rel = 1
        self.work_fluid.T1 = self.T_in_stag
        self.work_fluid.alpha = self.alpha
        while self.dk_rel >= 0.01:
            self.c_p_gas_old = self.work_fluid.c_p_av_int
            self.work_fluid.T2 = self.T_in_stag - self.L_turb / self.work_fluid.c_p_av_int
            logger.debug('comp_turb_output T_out_stag = %s' % self.work_fluid.T2)
            logger.debug('comp_turb_output c_p_gas_av = %s' % self.work_fluid.c_p_av_int)
            self.k_gas_new = self.work_fluid.k_av_int
            self.dk_rel = max(abs(self.k_gas_new - self.k_gas) / self.k_gas)
            self.k_gas_old = self.k_gas
            self.k_gas = self.k_gas_new
            logger.debug('comp_turb_output k = %s' % self.k_gas)
        self.pi_turb_stag = np.array([self.pi_turb_stag_init])
        self.dpi_rel = 1
        self.p_out_stag = None
        self.eta_turb_stag = None
        while self.dpi_rel >= 0.01:
            self.eta_turb_stag = func.eta_turb_stag(self.pi_turb_stag, self.k_gas, self.eta_turb_stag_p)
            self.p_out_stag = self.p_in_stag * (1 -
                                      self.L_turb / (self.work_fluid.c_p_av_int * self.eta_turb_stag * self.T_in_stag)) ** (
                                     self.k_gas / (self.k_gas - 1))
            logger.debug('comp_turb_output p_out_stag = %s' % self.p_out_stag)
            self.pi_turb_stag_new = self.p_in_stag / self.p_out_stag
            self.dpi_rel = max(abs(self.pi_turb_stag - self.pi_turb_stag_new) / self.pi_turb_stag)
            self.pi_turb_stag_old = self.pi_turb_stag
            self.pi_turb_stag = self.pi_turb_stag_new
        self.H_turb_stag = self.work_fluid.c_p_av_int * self.T_in_stag * (1 - self.pi_turb_stag ** ((1 - self.k_gas) / self.k_gas))
        return {'L_turb': self.L_turb, 'k': self.k_gas, 'c_p': self.work_fluid.c_p_av_int, 'T_out_stag': self.work_fluid.T2,
                'p_out_stag': self.p_out_stag, 'pi_turb_stag': self.pi_turb_stag, 'H_turb_stag': self.H_turb_stag,
                'eta_turb_stag': self.eta_turb_stag, 'work_fluid': self.work_fluid}

    @property
    def input(self) -> CompressorTurbineInput:
        return self._input

    @input.setter
    def input(self, value: CompressorTurbineInput):
        self._input = value

    @property
    def output(self):
        return self._output

    def compute_output(self):
        output = self._compute_comp_turb_output(self._input.work_fluid, self._input.pi_turb_stag_init, self._input.T_in_stag,
                                          self._input.p_in_stag, self._input.L_comp, self._input.g_gas,
                                          self._input.eta_turb_stag_p, self._input.eta_m, self._input.alpha)
        self._output.c_p_gas = output['c_p']
        self._output.H_turb_stag = output['H_turb_stag']
        self._output.k_gas = output['k']
        self._output.L_turb = output['L_turb']
        self._output.p_out_stag = output['p_out_stag']
        self._output.pi_turb_stag = output['pi_turb_stag']
        self._output.T_out_stag = output['T_out_stag']
        self._output.eta_turb_stag = output['eta_turb_stag']
        self._output.work_fluid = output['work_fluid']
        del(output['work_fluid'])
        self._output.output_frame = pandas.DataFrame.from_dict(output)


def compute_chamber_output(T_in_stag, T_out_stag, p_in_stag, sigma_comb, eta_comb, Q_n, l0, work_fluid_in: Air,
                           work_fluid_out: KeroseneCombustionProducts, g_outflow, g_cooling, g_return):
    work_fluid_out.__init__()
    work_fluid_in.__init__()
    work_fluid_in.T = T_in_stag
    work_fluid_out.T = T_out_stag
    work_fluid_out.alpha = 1
    work_fluid_out_clean = KeroseneCombustionProducts()
    work_fluid_out_clean.alpha = 1
    work_fluid_out_clean.T = 288
    g_fuel = (work_fluid_out.c_p_av * T_out_stag - work_fluid_in.c_p_av * T_in_stag) / \
             (Q_n * eta_comb - work_fluid_out.c_p_av * T_out_stag + work_fluid_out_clean.c_p_av * 288)
    alpha = 1 / (g_fuel * l0)
    g_gas = (1 + g_fuel) * (1 - g_outflow - g_cooling) + g_return
    p_out_stag = p_in_stag * sigma_comb
    return {'g_fuel': g_fuel, 'g_gas': g_gas, 'alpha': alpha, 'p_out_stag': p_out_stag}


class CombustionChamber:
    def __init__(self):
        self._input = CombustionChamberInput()
        self._output = CombustionChamberOutput()

    def _compute_chamber_output(self, T_in_stag, T_out_stag, p_in_stag, sigma_comb, eta_comb, Q_n, l0, work_fluid_in: Air,
                                work_fluid_out: KeroseneCombustionProducts, g_outflow, g_cooling, g_return):
        self.work_fluid_out = work_fluid_out
        self.work_fluid_out.__init__()
        self.work_fluid_in = work_fluid_in
        self.work_fluid_in.__init__()
        self.T_in_stag = T_in_stag
        self.T_out_stag = T_out_stag
        self.p_in_stag = p_in_stag
        self.sigma_comb = sigma_comb
        self.eta_comb = eta_comb
        self.Q_n = Q_n
        self.l0 = l0
        self.g_outflow = g_outflow
        self.g_cooling = g_cooling
        self.g_return = g_return
        self.work_fluid_in.T = T_in_stag
        self.work_fluid_out.T = T_out_stag
        self.work_fluid_out.alpha = 1
        self.work_fluid_out_clean = KeroseneCombustionProducts()
        self.work_fluid_out_clean.alpha = 1
        self.work_fluid_out_clean.T = 288
        self.g_fuel = (self.work_fluid_out.c_p_av * self.T_out_stag - self.work_fluid_in.c_p_av * self.T_in_stag) / \
                 (self.Q_n * self.eta_comb - self.work_fluid_out.c_p_av * self.T_out_stag + self.work_fluid_out_clean.c_p_av * 288)
        self.alpha = 1 / (self.g_fuel * self.l0)
        self.g_gas = (1 + self.g_fuel) * (1 - self.g_outflow - self.g_cooling) + self.g_return
        self.p_out_stag = self.p_in_stag * self.sigma_comb
        return {'g_fuel': self.g_fuel, 'g_gas': self.g_gas, 'alpha': self.alpha, 'p_out_stag': self.p_out_stag}

    @property
    def input(self) -> CombustionChamberInput:
        return self._input

    @input.setter
    def input(self, value: CombustionChamberInput):
        self._input = value

    @property
    def output(self) -> CombustionChamberOutput:
        return self._output

    def compute_output(self):
        output = self._compute_chamber_output(self._input.T_in_stag, self._input.T_out_stag, self._input.p_in_stag,
                                        self._input.sigma_comb, self._input.eta_comb, self._input.Q_n, self._input.l0,
                                        self._input.work_fluid_in, self._input.work_fluid_out, self._input.g_outflow,
                                        self._input.g_cooling, self._input.g_return)
        self._output.alpha = output['alpha']
        self._output.g_fuel = output['g_fuel']
        self._output.g_gas = output['g_gas']
        self._output.p_out_stag = output['p_out_stag']
        self._output.output_frame = pandas.DataFrame.from_dict(output)


def compute_inlet_output(work_fluid: Air, T_in, p_in, M_v, v, sigma_in):
    work_fluid.__init__()
    T_out_stag = None
    p_out_stag = None
    pi_v = None
    if (M_v == 0) or (v == 0):
        work_fluid.T = T_in
        T_out_stag = T_in
        p_out_stag = p_in * sigma_in
        pi_v = p_out_stag / p_in
    else:
        if v is None:
            work_fluid.T1 = T_in
            k_air = work_fluid.k_av_int
            dk_rel = 1
            while dk_rel >= 0.01:
                T_out_stag = T_in * gd.GasDynamicFunctions.tau_M(M_v, k_air)
                work_fluid.T2 = T_out_stag
                k_air_new = work_fluid.k_av_int
                dk_rel = abs(k_air_new - k_air) / k_air
                k_air = k_air_new
            p_out_stag = p_in * gd.GasDynamicFunctions.pi_M(M_v, k_air) * sigma_in
            pi_v = p_out_stag / p_in
        elif M_v is None or (M_v is not None and v is not None):
            work_fluid.T1 = T_in
            k_air = work_fluid.k_av_int
            dk_rel = 1
            while dk_rel >= 0.01:
                T_out_stag = T_in + v**2 / (2 * work_fluid.c_p_av_int)
                work_fluid.T2 = T_out_stag
                k_air_new = work_fluid.k_av_int
                dk_rel = abs(k_air_new - k_air) / k_air
                k_air = k_air_new
            p_out_stag = p_in * (1 + v**2 / (2 * work_fluid.c_p_av_int * T_in))**\
                                (work_fluid.k_av_int / (work_fluid.k_av_int - 1)) * sigma_in
            pi_v = p_out_stag / p_in
    return {'T_out_stag': T_out_stag, 'p_out_stag': p_out_stag, 'pi_v': pi_v}


class Inlet:
    def __init__(self):
        self._input = InletInput()
        self._output = InletOutput()

    def _compute_inlet_output(self, work_fluid: Air, T_in, p_in, M_v, v, sigma_in):
        self.work_fluid = work_fluid
        self.work_fluid.__init__()
        self.T_in = T_in
        self.M_v = M_v
        self.p_in = p_in
        self.v = v
        self.sigma_in = sigma_in
        self.T_out_stag = None
        self.p_out_stag = None
        self.pi_v = None
        if (self.M_v == 0) or (self.v == 0):
            self.work_fluid.T = self.T_in
            self.T_out_stag = self.T_in
            self.p_out_stag = self.p_in * self.sigma_in
            self.pi_v = self.p_out_stag / self.p_in
        else:
            if self.v is None:
                self.work_fluid.T1 = self.T_in
                self.k_air = self.work_fluid.k_av_int
                self.dk_rel = 1
                while self.dk_rel >= 0.01:
                    self.T_out_stag = self.T_in * gd.GasDynamicFunctions.tau_M(self.M_v, self.k_air)
                    self.work_fluid.T2 = self.T_out_stag
                    self.k_air_new = self.work_fluid.k_av_int
                    self.dk_rel = abs(self.k_air_new - self.k_air) / self.k_air
                    self.k_air = self.k_air_new
                self.p_out_stag = self.p_in * gd.GasDynamicFunctions.pi_M(self.M_v, self.k_air) * self.sigma_in
                self.pi_v = self.p_out_stag / self.p_in
            elif self.M_v is None or (self.M_v is not None and self.v is not None):
                self.work_fluid.T1 = self.T_in
                self.k_air = self.work_fluid.k_av_int
                self.dk_rel = 1
                while self.dk_rel >= 0.01:
                    self.T_out_stag = self.T_in + self.v ** 2 / (2 * self.work_fluid.c_p_av_int)
                    self.work_fluid.T2 = self.T_out_stag
                    self.k_air_new = self.work_fluid.k_av_int
                    self.dk_rel = abs(self.k_air_new - self.k_air) / self.k_air
                    self.k_air = self.k_air_new
                self.p_out_stag = self.p_in * (1 + self.v ** 2 / (2 * self.work_fluid.c_p_av_int * self.T_in)) ** \
                                    (self.work_fluid.k_av_int / (self.work_fluid.k_av_int - 1)) * self.sigma_in
                self.pi_v = self.p_out_stag / self.p_in
        return {'T_out_stag': self.T_out_stag, 'p_out_stag': self.p_out_stag, 'pi_v': self.pi_v}

    @property
    def input(self) -> InletInput:
        return self._input

    @input.setter
    def input(self, value: InletInput):
        self._input = value

    @property
    def output(self) -> InletOutput:
        return self._output

    def compute_output(self):
        output = self._compute_inlet_output(self._input.work_fluid, self._input.T_in, self._input.p_in, self._input.M_v,
                                      self._input.v, self._input.sigma_in)
        self._output.p_out_stag = output['p_out_stag']
        self._output.pi_v = output['pi_v']
        self._output.T_out_stag = output['T_out_stag']


def compute_outlet_output(work_fluid: KeroseneCombustionProducts, T_in_stag, lam_out, p_out, sigma_out, alpha):
    work_fluid.__init__()
    work_fluid.alpha = alpha
    work_fluid.T = T_in_stag
    gd_par = gd.GasDynamicsParameters(k=work_fluid.k, R=work_fluid.R, lam=lam_out, p=p_out, T_stag=T_in_stag)
    p_in_stag = gd_par.p_stag / sigma_out
    return {'T_out_stag': T_in_stag, 'p_in_stag': p_in_stag, 'T_out': gd_par.T, 'p_out_stag': gd_par.p_stag,
            'c_out': gd_par.c, 'k': work_fluid.k, 'c_p': work_fluid.c_p}


class Outlet:
    def __init__(self):
        self._input = OutletInput()
        self._output = OutletOutput()

    def _compute_outlet_output(self, work_fluid: KeroseneCombustionProducts, T_in_stag, lam_out, p_out, sigma_out, alpha):
        self.work_fluid = work_fluid
        self.work_fluid.__init__()
        self.T_in_stag = T_in_stag
        self.lam_out = lam_out
        self.p_out = p_out
        self.sigma_out = sigma_out
        self.alpha = alpha
        self.work_fluid.alpha = alpha
        self.work_fluid.T = T_in_stag
        self.gd_par = gd.GasDynamicsParameters(k=self.work_fluid.k, R=self.work_fluid.R, lam=self.lam_out,
                                          p=self.p_out, T_stag=self.T_in_stag)
        self.p_in_stag = self.gd_par.p_stag / self.sigma_out
        return {'T_out_stag': self.T_in_stag, 'p_in_stag': self.p_in_stag, 'T_out': self.gd_par.T,
                'p_out_stag': self.gd_par.p_stag, 'c_out': self.gd_par.c, 'k': self.work_fluid.k,
                'c_p': self.work_fluid.c_p}

    @property
    def input(self) -> OutletInput:
        return self._input

    @input.setter
    def input(self, value: OutletInput):
        self._input = value

    @property
    def output(self) -> OutletOutput:
        return self._output

    def compute_output(self):
        output = self._compute_outlet_output(self._input.work_fluid, self._input.T_in_stag, self._input.lam_out,
                                       self._input.p_out, self._input.sigma_out, self._input.alpha)
        self._output.c_out = output['c_out']
        self._output.p_in_stag = output['p_in_stag']
        self._output.c_p_gas = output['c_p']
        self._output.k_gas = output['k']
        self._output.T_out_stag = output['T_out_stag']
        self._output.T_out = output['T_out']
        self._output.p_out_stag = output['p_out_stag']
        self._output.output_frame = pandas.DataFrame.from_dict(output)


if __name__ == '__main__':
    compressor = Compressor()
    compressor.input.eta_comp_stag_p = 0.9
    compressor.input.pi_comp_stag = np.arange(5, 20, 1)
    compressor.input.work_fluid = Air()
    compressor.input.T_in_stag = 300
    compressor.input.p_in_stag = 1e5
    compressor.compute_output()
    print(compressor.output.output_frame)
    print()

    power_turbine = PowerTurbine()
    power_turbine.input.alpha = np.array([2, 2.2, 1.5])
    power_turbine.input.c_out = 200
    power_turbine.input.eta_turb_stag_p = 0.9
    power_turbine.input.T_out_stag_init = 900
    power_turbine.input.work_fluid = KeroseneCombustionProducts()
    power_turbine.input.p_in_stag = np.array([500000, 600000, 800000])
    power_turbine.input.T_in_stag = np.array([1400, 1450, 1500])
    power_turbine.input.p_out_stag = 1e5
    power_turbine.compute_output()
    print(power_turbine.output.output_frame)
    print()

    comp_turbine = CompressorTurbine()
    comp_turbine.input.alpha = np.array([2, 2.2, 1.5])
    comp_turbine.input.c_out = 220
    comp_turbine.input.eta_m = 0.99
    comp_turbine.input.eta_turb_stag_p = 0.91
    comp_turbine.input.pi_turb_stag_init = 3
    comp_turbine.input.g_gas = 0.95
    comp_turbine.input.work_fluid = KeroseneCombustionProducts()
    comp_turbine.input.L_comp = 4e5
    comp_turbine.input.p_in_stag = np.array([6e5, 7e5, 8e5])
    comp_turbine.input.T_in_stag = np.array([1300, 1400, 1700])
    comp_turbine.compute_output()
    print(comp_turbine.output.output_frame)
