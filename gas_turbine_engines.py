import cycle_calculation.gas_turbine_units as gt_units
import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
import pickle as pk
import os


logger = gt_units.func.create_logger(loggerlevel=logging.INFO, filename=gt_units.log_file_name, filemode='w',
                                     name=__name__, add_console_handler=False)


class TV7_117:

    _output_fields = ['_N_e_specific', '_G_air', '_C_e', '_eta_e', '_pi_comp_stag_L', '_pi_comp_stag_eta',
                      '_pi_comp_stag_opt', '_compressor', '_combustion_chamber', '_compressor_turbine', '_inlet',
                      '_power_turbine', '_outlet']

    def __init__(self):
        self.T_comb_stag = None
        self.T_a = None
        self.p_a = None
        self.M_v = None
        self.pi_comp_stag = None
        self.sigma_in = None
        self.sigma_out = None
        self.sigma_comb = None
        self.N_e = None
        self.eta_comp_stag_p = None
        self.eta_comp_turb_stag_p = None
        self.eta_comb = None
        self.eta_power_turb_stag_p = None
        self.eta_m = None
        self.eta_r = None
        self.Q_n = None
        self.l0 = None
        self.g_outflow = None
        self.g_cooling = None
        self.g_return = None
        self.lam_out = None
        self.comp_work_fluid = None
        self.turb_work_fluid = None
        self._N_e_specific = None
        self._G_air = None
        self._C_e = None
        self._eta_e = None
        self._pi_comp_stag_L = None
        self._pi_comp_stag_eta = None
        self._pi_comp_stag_opt = None
        self._inlet = gt_units.Inlet()
        self._compressor = gt_units.Compressor()
        self._combustion_chamber = gt_units.CombustionChamber()
        self._compressor_turbine = gt_units.CompressorTurbine()
        self._power_turbine = gt_units.PowerTurbine()
        self._outlet = gt_units.Outlet()

    @property
    def outlet(self) -> gt_units.Outlet:
        return self._outlet

    @property
    def N_e_specific(self):
        return self._N_e_specific

    @property
    def G_air(self):
        return self._G_air

    @property
    def inlet(self) -> gt_units.Inlet:
        return self._inlet

    @property
    def compressor(self) -> gt_units.Compressor:
        return self._compressor

    @property
    def combustion_chamber(self) -> gt_units.CombustionChamber:
        return self._combustion_chamber

    @property
    def compressor_turbine(self) -> gt_units.CompressorTurbine:
        return self._compressor_turbine

    @property
    def power_turbine(self) -> gt_units.PowerTurbine:
        return self._power_turbine

    @property
    def C_e(self):
        return self._C_e

    @property
    def eta_e(self):
        return self._eta_e

    @property
    def pi_comp_stag_L(self):
        return self._pi_comp_stag_L

    @property
    def pi_comp_stag_eta(self):
        return self._pi_comp_stag_eta

    @property
    def pi_comp_stag_opt(self):
        return self._pi_comp_stag_opt

    def __str__(self):
        return 'TV7_117 object'

    def _compute_inlet(self):
        self._inlet.input.work_fluid = self.comp_work_fluid
        self._inlet.input.M_v = self.M_v
        self._inlet.input.p_in = self.p_a
        self._inlet.input.sigma_in = self.sigma_in
        self._inlet.input.T_in = self.T_a
        logger.info('%s _compute_inlet' % str(self))
        self._inlet.compute_output()

    def _compute_compressor(self):
        self._compressor.input.eta_comp_stag_p = self.eta_comp_stag_p
        self._compressor.input.p_in_stag = self._inlet.output.p_out_stag
        self._compressor.input.pi_comp_stag = self.pi_comp_stag
        self._compressor.input.T_in_stag = self._inlet.output.T_out_stag
        self._compressor.input.work_fluid = self.comp_work_fluid
        logger.info('%s _compute_compressor' % str(self))
        self._compressor.compute_output()

    def _compute_combustion_chamber(self):
        self._combustion_chamber.input.T_out_stag = self.T_comb_stag
        self._combustion_chamber.input.eta_comb = self.eta_comb
        self._combustion_chamber.input.g_cooling = self.g_cooling
        self._combustion_chamber.input.g_outflow = self.g_outflow
        self._combustion_chamber.input.g_return = self.g_return
        self._combustion_chamber.input.l0 = self.l0
        self._combustion_chamber.input.p_in_stag = self._compressor.output.p_out_stag
        self._combustion_chamber.input.Q_n = self.Q_n
        self._combustion_chamber.input.sigma_comb = self.sigma_comb
        self._combustion_chamber.input.T_in_stag = self._compressor.output.T_out_stag
        self._combustion_chamber.input.work_fluid_in = self.comp_work_fluid
        self._combustion_chamber.input.work_fluid_out = self.turb_work_fluid
        logger.info('%s _compute_combustion_chamber' % str(self))
        self._combustion_chamber.compute_output()

    def _compute_compressor_turbine(self):
        self._compressor_turbine.input.T_in_stag = self.T_comb_stag
        self._compressor_turbine.input.p_in_stag = self._combustion_chamber.output.p_out_stag
        self._compressor_turbine.input.alpha = self._combustion_chamber.output.alpha
        self._compressor_turbine.input.eta_m = self.eta_m
        self._compressor_turbine.input.eta_turb_stag_p = self.eta_comp_turb_stag_p
        self._compressor_turbine.input.g_gas = self._combustion_chamber.output.g_gas
        self._compressor_turbine.input.L_comp = self._compressor.output.L_comp
        self._compressor_turbine.input.pi_turb_stag_init = 3
        self._compressor_turbine.input.work_fluid = self.turb_work_fluid
        logger.info('%s _compute_compressor_turbine' % str(self))
        self._compressor_turbine.compute_output()

    def _compute_power_turbine(self):
        self._power_turbine.input.alpha = self._combustion_chamber.output.alpha
        self._power_turbine.input.work_fluid = self.turb_work_fluid
        self._power_turbine.input.eta_turb_stag_p = self.eta_power_turb_stag_p
        self._power_turbine.input.p_in_stag = self._compressor_turbine.output.p_out_stag
        self._power_turbine.input.T_in_stag = self._compressor_turbine.output.T_out_stag
        self._power_turbine.input.T_out_stag_init = 600
        logger.info('%s _compute_power_turbine' % str(self))
        self._power_turbine.compute_output()

    def _compute_outlet(self):
        self._outlet.input.alpha = self._combustion_chamber.output.alpha
        self._outlet.input.lam_out = self.lam_out
        self._outlet.input.p_out = self.p_a
        self._outlet.input.sigma_out = self.sigma_out
        self._outlet.input.work_fluid = self.turb_work_fluid
        logger.info('%s _compute_outlet' % str(self))
        self._outlet.compute_output()

    def _compute_power_turbine_outlet_complex(self):
        self._power_turbine.input.p_out_stag = np.array([self.p_a])
        delta = 1
        while delta >= 0.01:
            self._compute_power_turbine()
            self._outlet.input.T_in_stag = self._power_turbine.output.T_out_stag
            self._compute_outlet()
            delta = max(abs(self._outlet.output.p_in_stag -
                        self._power_turbine.input.p_out_stag) / self._power_turbine.input.p_out_stag)
            self._power_turbine.input.p_out_stag = self._outlet.output.p_in_stag

    def _compute_optimal_pi(self):
        set_atr = object.__setattr__
        if len(self.pi_comp_stag) > 1:
            logger.info('%s _compute_optimal_pi' % str(self))
            G_air_interp = interp1d(self.pi_comp_stag, self._G_air)

            def G_air(pi_comp_stag):
                return G_air_interp(pi_comp_stag)

            C_e_interp = interp1d(self.pi_comp_stag, self._C_e)

            def C_e(pi_comp_stag):
                return C_e_interp(pi_comp_stag)

            pi_comp_stag_L = minimize_scalar(G_air,
                                             bounds=(self.pi_comp_stag[0], self.pi_comp_stag[len(self.pi_comp_stag) - 1]),
                                             method='bounded')['x']
            pi_comp_stag_eta = minimize_scalar(C_e,
                                               bounds=(self.pi_comp_stag[0], self.pi_comp_stag[len(self.pi_comp_stag) - 1]),
                                               method='bounded')['x']
            pi_comp_stag_opt = 0.35 * (pi_comp_stag_eta - pi_comp_stag_L) + pi_comp_stag_L
        else:
            pi_comp_stag_L = self.pi_comp_stag[0]
            pi_comp_stag_eta = self.pi_comp_stag[0]
            pi_comp_stag_opt = self.pi_comp_stag[0]
        set_atr(self, '_pi_comp_stag_L', np.array([pi_comp_stag_L]))
        set_atr(self, '_pi_comp_stag_eta', np.array([pi_comp_stag_eta]))
        set_atr(self, '_pi_comp_stag_opt', np.array([pi_comp_stag_opt]))

    def _compute_all_output(self):
        set_atr = object.__setattr__
        self._compute_inlet()
        self._compute_compressor()
        self._compute_combustion_chamber()
        self._compute_compressor_turbine()
        self._compute_power_turbine_outlet_complex()
        N_e_specific_value = (self._power_turbine.output.L_turb *
                              self._combustion_chamber.output.g_gas * self.eta_m * self.eta_r)
        G_air_value = self.N_e / N_e_specific_value
        C_e = self._combustion_chamber.output.g_fuel * (1 - self.g_cooling -
                                                        self.g_outflow) * 3600 / N_e_specific_value
        eta_e = 3600 / (C_e * self.Q_n)
        set_atr(self, '_N_e_specific', N_e_specific_value)
        set_atr(self, '_G_air', G_air_value)
        set_atr(self, '_C_e', C_e)
        set_atr(self, '_eta_e', eta_e)
        self._compute_optimal_pi()

    def _none_value_count(self):
        k = 0
        for item in self.__dict__:
            if (self.__dict__[item] is None) and (item not in self._output_fields):
                k += 1
        return k

    def __setattr__(self, key, value):
        set_atr = object.__setattr__
        get_atr = object.__getattribute__
        if (key not in self._output_fields) and (self._none_value_count() != 0):
            set_atr(self, key, value)
            if self._none_value_count() == 0:
                self._compute_all_output()
        elif (key not in self._output_fields) and (self._none_value_count() == 0):
            set_atr(self, key, value)
            try:
                self._compute_all_output()
            except AttributeError:
                pass
        elif key in self._output_fields:
            try:
                get_atr(self, key)
            except AttributeError:
                set_atr(self, key, value)

    def plot(self, y_arr: np.ndarray, ylabel):
        plt.plot(self.pi_comp_stag, y_arr, linewidth=2)
        plt.grid()
        plt.xlabel(r'$\pi_k$', fontsize=18)
        plt.ylabel(ylabel, fontsize=18)
        plt.show()

    def save_output(self, filename='cycle_calculation_results'):
        file = open(os.path.join(os.path.dirname(__file__), filename), 'wb')
        pk.dump(self, file)
        file.close()

