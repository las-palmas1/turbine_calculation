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
    while dk_rel >= 0.01:
        k += 1
        eta_comp_stag = func.eta_comp_stag(pi_comp_stag, k_air_ret, eta_comp_stag_p)
        work_fluid.T2 = T_in_stag * (1 + (pi_comp_stag ** ((k_air_ret - 1) / k_air_ret) - 1) / eta_comp_stag)
        k_air_ret_new = work_fluid.k_av_int
        dk_rel = max(abs(k_air_ret - k_air_ret_new) / k_air_ret)
        k_air_ret = k_air_ret_new
    L_comp = work_fluid.c_p_av_int * (work_fluid.T2 - work_fluid.T1)
    return {'k_air': k_air_ret, 'c_p_air': work_fluid.c_p_av_int, 'eta_comp_stag': eta_comp_stag,
            'T_out_stag': work_fluid.T2, 'L_comp': L_comp, 'work_fluid': work_fluid}


class Compressor:
    def __init__(self):
        self._input = CompressorInput()
        self._output = CompressorOutput()

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
        output = compute_comp_output(self._input.work_fluid, self._input.pi_comp_stag,
                                     self._input.eta_comp_stag_p, self._input.T_in_stag)
        self._output.k_air = output['k_air']
        self._output.p_out_stag = self._input.p_in_stag * self._input.pi_comp_stag
        self._output.c_p_air = output['c_p_air']
        self._output.eta_comp_stag = output['eta_comp_stag']
        self._output.T_out_stag = output['T_out_stag']
        self._output.L_comp = output['L_comp']
        self._output.work_fluid = output['work_fluid']
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
        work_fluid.T2 = T_in_stag * (1 - (1 - pi_turb_stag **
                                          ((1 - work_fluid.k_av_int) / work_fluid.k_av_int)) * eta_turb_stag)
        dT_rel = max(abs(work_fluid.T2 - work_fluid.T2) / work_fluid.T2)
    L_turb = work_fluid.c_p_av_int * (T_in_stag - work_fluid.T2)
    H_turb_stag = work_fluid.c_p_av_int * \
                  (1 - pi_turb_stag**((1 - work_fluid.k_av_int) / work_fluid.k_av_int)) * T_in_stag
    return {'k_gas': work_fluid.k_av_int, 'c_p_gas': work_fluid.c_p_av_int, 'pi_turb_stag': pi_turb_stag,
            'T_out_stag': work_fluid.T2, 'L_turb': L_turb, 'H_turb_stag': H_turb_stag,
            'eta_turb_stag': eta_turb_stag, 'work_fluid': work_fluid}


class PowerTurbine:
    def __init__(self):
        self._input = PowerTurbineInput()
        self._output = PowerTurbineOutput()

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
        output = compute_power_turb_output(self._input.work_fluid, self._input.T_out_stag_init, self._input.T_in_stag,
                                           self._input.p_in_stag, self.input.p_out_stag, self._input.eta_turb_stag_p,
                                           self._input.alpha)
        self._output.k_gas = output['k_gas']
        self._output.c_p_gas = output['c_p_gas']
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
        k_gas = k_gas_new
        logger.debug('comp_turb_output k_gas = %s' % k_gas)
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
    return {'L_turb': L_turb, 'k_gas': k_gas, 'c_p_gas': work_fluid.c_p_av_int, 'T_out_stag': work_fluid.T2,
            'p_out_stag': p_out_stag, 'pi_turb_stag': pi_turb_stag, 'H_turb_stag': H_turb_stag,
            'eta_turb_stag': eta_turb_stag, 'work_fluid': work_fluid}


class CompressorTurbine:
    def __init__(self):
        self._input = CompressorTurbineInput()
        self._output = CompressorTurbineOutput()

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
        output = compute_comp_turb_output(self._input.work_fluid, self._input.pi_turb_stag_init, self._input.T_in_stag,
                                          self._input.p_in_stag, self._input.L_comp, self._input.g_gas,
                                          self._input.eta_turb_stag_p, self._input.eta_m, self._input.alpha)
        self._output.c_p_gas = output['c_p_gas']
        self._output.H_turb_stag = output['H_turb_stag']
        self._output.k_gas = output['k_gas']
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
    g_fuel = (work_fluid_out.c_p * T_out_stag - work_fluid_in.c_p * T_in_stag) / \
             (Q_n * eta_comb - work_fluid_out.c_p * T_out_stag + work_fluid_out_clean.c_p * 288)
    alpha = 1 / (g_fuel * l0)
    g_gas = (1 + g_fuel) * (1 - g_outflow - g_cooling) + g_return
    p_out_stag = p_in_stag * sigma_comb
    return {'g_fuel': g_fuel, 'g_gas': g_gas, 'alpha': alpha, 'p_out_stag': p_out_stag}


class CombustionChamber:
    def __init__(self):
        self._input = CombustionChamberInput()
        self._output = CombustionChamberOutput()

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
        output = compute_chamber_output(self._input.T_in_stag, self._input.T_out_stag, self._input.p_in_stag,
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
        output = compute_inlet_output(self._input.work_fluid, self._input.T_in, self._input.p_in, self._input.M_v,
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
            'c_out': gd_par.c, 'k_gas': work_fluid.k, 'c_p_gas': work_fluid.c_p}


class Outlet:
    def __init__(self):
        self._input = OutletInput()
        self._output = OutletOutput()

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
        output = compute_outlet_output(self._input.work_fluid, self._input.T_in_stag, self._input.lam_out,
                                       self._input.p_out, self._input.sigma_out, self._input.alpha)
        self._output.c_out = output['c_out']
        self._output.p_in_stag = output['p_in_stag']
        self._output.c_p_gas = output['c_p_gas']
        self._output.k_gas = output['k_gas']
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
