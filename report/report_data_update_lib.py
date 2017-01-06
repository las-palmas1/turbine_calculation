from cycle_calculation.gas_turbine_engines import TV7_117
from average_streamline_calculation.turbine import Turbine
import pickle as pc
import config
import strength_calculation.config as strength_config
from strength_calculation.post_processing import two_calculations, material, p_m
from parts_parameters.disks import FirstStageDisk
import numpy as np


def get_TV7_117_object(filename) -> TV7_117:
    file = open(filename, 'rb')
    result = pc.load(file)
    file.close()
    return result


def get_turbine_object(filename) -> Turbine:
    file = open(filename, 'rb')
    result = pc.load(file)
    file.close()
    return result


class NumInserter:
    def __init__(self, input_path, output_path, parameter_dict):
        self.input_path = input_path
        self.output_path = output_path
        self.parameter_dict = parameter_dict

    def process_document(self):
        content = self._get_file_content()
        modified_content = self._modify_fle_content(content)
        self._write_modified_content(modified_content)

    def _write_modified_content(self, content: str):
        with open(self.output_path, 'w', encoding='utf-8') as output_file:
            output_file.write(content)

    def _modify_fle_content(self, content: str):
        result = content
        for key in self.parameter_dict:
            result = result.replace(r'<%s>' % key, str(self.parameter_dict[key]))

        return result

    def _get_file_content(self):
        file = open(self.input_path, encoding='utf-8')
        content = file.read()
        file.close()

        return content


def get_parameter_dict(engine: TV7_117, turbine: Turbine):
    def get_first_calculation_rows():
        template = r'%.0f & %.1f & %.1f & %.1f \\ \hline '
        res = ''
        for n in range(two_calculations.part_number):
            res += template % (round(n, 0), round(two_calculations.first_calc.r_arr[n] * 1e3, 1),
                               round(two_calculations.first_calc.sigma_r_arr[n][0] / 1e6, 1),
                               round(two_calculations.first_calc.sigma_t_arr[n][0] / 1e6, 1))
            res += template % (round(n, 0), round(two_calculations.first_calc.r_arr[n + 1] * 1e3, 1),
                               round(two_calculations.first_calc.sigma_r_arr[n][1] / 1e6, 1),
                               round(two_calculations.first_calc.sigma_t_arr[n][1] / 1e6, 1))
        return res

    def get_second_calculation_rows():
        template = r'%.0f & %f & %f & %f \\ \hline '
        res = ''
        for n in range(two_calculations.part_number):
            res += template % (round(n, 0), round(two_calculations.second_calc.r_arr[n] * 1e3, 1),
                               round(two_calculations.second_calc.sigma_r_arr[n][0] / 1e6, 1),
                               round(two_calculations.second_calc.sigma_t_arr[n][0] / 1e6, 1))
            res += template % (round(n, 0), round(two_calculations.second_calc.r_arr[n + 1] * 1e3, 1),
                               round(two_calculations.second_calc.sigma_r_arr[n][1] / 1e6, 1),
                               round(two_calculations.second_calc.sigma_t_arr[n][1] / 1e6, 1))
        return res

    def get_recalculated_stress_rows():
        template = r'%.0f & %f & %f & %f \\ \hline '
        res = ''
        for n in range(two_calculations.part_number):
            res += template % (round(n, 0), round(two_calculations.r_arr[n] * 1e3, 1),
                               round(two_calculations.sigma_r_arr_real[n][0] / 1e6, 1),
                               round(two_calculations.sigma_t_arr_real[n][0] / 1e6, 1))
            res += template % (round(n, 0), round(two_calculations.r_arr[n + 1] * 1e3, 1),
                               round(two_calculations.sigma_r_arr_real[n][1] / 1e6, 1),
                               round(two_calculations.sigma_r_arr_real[n][1] / 1e6, 1))
        return res

    def get_average_stress_rows():
        template = r'%.0f & %f & %f & %f & %f \\ \hline '
        res = ''
        for n in range(two_calculations.part_number):
            res += template % (round(n, 0), round(two_calculations.r_arr[n] * 1e3, 1),
                               round(two_calculations.sigma_r_arr_av[n] / 1e6, 1),
                               round(two_calculations.sigma_t_arr_av[n] / 1e6, 1),
                               round(two_calculations.sigma_eq_arr[n] / 1e6, 1))
        return res

    result = {
        # Исходные данные
        'EtaCompPol': config.eta_comp_stag_p,
        'EtaBurn': config.eta_comb,
        'EtaCompTurbPol': config.eta_comp_turb_stag_p,
        'EtaFreeTurbPol': config.eta_power_turb_stag_p,
        'LambdaOut': config.lam_out,
        'CyclePower': round(config.N_e / 1e3, 0),
        'GasTemp': config.T_comb_stag,
        'InletPipeSigma': config.sigma_in,
        'OutletPipeSigma': config.sigma_out,
        'EtaMech': config.eta_m,
        'EtaR': config.eta_r,
        'LossMassRateRel': config.g_outflow,
        'CoolMassRateRel': config.g_cooling,
        'ReturnMassRateRel': config.g_return,

        # Расчет цикла
        'CyclePlotGAir': r'\begin{figure}[hbtp] '
                         r'\centering '
                         r'\includegraphics[scale=0.5] {../../plots/cycle_G_air.png} '
                         r'\caption{Зависимость расхода воздуха от степени повышения давления} '
                         r'\end{figure}',
        'CyclePlotCe': r'\begin{figure}[hbtp] '
                         r'\centering '
                         r'\includegraphics[scale=0.5] {../../plots/cycle_C_e.png} '
                         r'\caption{Зависимость удельного расхода топлива от степени повышения давления} '
                         r'\end{figure}',
        'CyclePlotEtaE': r'\begin{figure}[hbtp] '
                       r'\centering '
                       r'\includegraphics[scale=0.5] {../../plots/cycle_eta_e.png} '
                       r'\caption{Зависимость КПД двигателя от степени повышения давления} '
                       r'\end{figure}',
        'PiComp': round(engine.pi_comp_stag_opt[0], 2),
        'AtmP': round(config.p_a / 1e6, 4),
        'AtmT': round(config.T_a, 0),
        'InletTubePOut': round(engine.inlet.output.p_out_stag / 1e6, 3),
        'CompPOut': round(engine.compressor.output.p_out_stag[0] / 1e6, 3),
        'CycleKAirShort': round(engine.compressor.k_air_old[0], 3),
        'EtaComp': round(engine.compressor.eta_comp_stag[0], 4),
        'CompTOut': round(engine.compressor.output.T_out_stag[0], 3),
        'CycleKAirLong': round(engine.compressor.k_air_ret[0], 3),
        'CycleKAirDelta': round(engine.compressor.dk_rel * 100, 2),
        'CycleAirGasConstant': engine.compressor.work_fluid.R,
        'CycleAirSpecificHeat': round(engine.compressor.work_fluid.c_p_av_int[0], 3),
        'CompSpecificLabour': round(engine.compressor.L_comp[0] / 1e6, 3),
        'BurnerTOut': round(engine.combustion_chamber.T_out_stag, 3),
        'ParameterDetermT': 290,
        'BurnerInletAirSpecificHeat': round(engine.combustion_chamber.work_fluid_in.c_p_av[0], 3),
        'BurnerOutletGasSpecificHeat': round(engine.combustion_chamber.work_fluid_out.c_p_av, 3),
        'ParameterDetermGasSpecificHeat': round(engine.combustion_chamber.work_fluid_out_clean.c_p_av, 3),
        'LowerQ': round(config.Q_n / 1e3, 0),
        'TheoryAirMass': config.l0,
        'FuelMassRateRel': round(engine.combustion_chamber.g_fuel[0], 4),
        'CycleBurnAlphaLong': round(engine.combustion_chamber.alpha[0], 3),
        'CompTurbineMassRateRel': round(engine.combustion_chamber.g_gas[0], 3),
        'CompTurbineSpecificLabour': round(engine.compressor_turbine.L_turb[0] / 1e6, 3),
        'BurnSigma': config.sigma_comb,
        'CompTurbinePIn': round(engine.compressor_turbine.p_in_stag[0] / 1e6, 3),
        'CycleCompTurbineKGasShort': round(engine.compressor_turbine.k_gas_old[0], 3),
        'CycleCompTurbineGasConstant': engine.compressor_turbine.work_fluid.R,
        'CompTurbineSpecificHeat': round(engine.compressor_turbine.c_p_gas_old[0], 3),
        'CompTurbineTOut': round(engine.compressor_turbine.output.T_out_stag[0], 3),
        'CycleCompTurbineKGasLong': round(engine.compressor_turbine.k_gas_new[0], 3),
        'CompTurbineKCalcError': round(engine.compressor_turbine.dk_rel * 100, 2),
        'CycleCompTurbinePiShort': round(engine.compressor_turbine.pi_turb_stag_old[0], 3),
        'EtaCompTurb': round(engine.compressor_turbine.eta_turb_stag[0], 3),
        'CompTurbinePOut': round(engine.compressor_turbine.p_out_stag[0] / 1e6, 3),
        'CycleCompTurbinePiLong': round(engine.compressor_turbine.pi_turb_stag[0], 3),
        'CompTurbinePiCalcError': round(engine.compressor_turbine.dpi_rel * 100, 2),
        'CycleFreeTurbineKGasShort': round(engine. power_turbine.k_gas_old[0], 3),
        'OutletTubePOut': round(engine.outlet.output.p_out_stag[0] / 1e6, 3),
        'FreeTurbinePOut': round(engine.outlet.output.p_in_stag[0] / 1e6, 3),
        'FreeTurbinePi': round(engine.power_turbine.pi_turb_stag[0], 3),
        'EtaFreeTurb': round(engine.power_turbine.eta_turb_stag[0], 4),
        'FreeTurbineTOut': round(engine.power_turbine.output.T_out_stag[0]),
        'CycleFreeTurbineKGasLong': round(engine.power_turbine.k_gas_new[0], 3),
        'FreeTurbineKCalcError': round(engine.power_turbine.dk_rel, 3),
        'CycleFreeTurbineGasGasConstant': engine.power_turbine.work_fluid.R,
        'FreeTurbineSpecificHeat': round(engine.power_turbine.work_fluid.c_p_av_int[0], 3),
        'FreeTurbineSpecificLabour': round(engine.power_turbine.L_turb[0] / 1e6, 3),
        'FreeTurbineMassRateRel': round(engine.combustion_chamber.g_gas[0], 3),
        'EngineSpecificPower': round(engine.N_e_specific[0] / 1e6, 4),
        'EngineFuelMassRateSpecific': round(engine.C_e[0] * 1e3, 4),
        'EngineEta': round(engine.eta_e[0], 4),
        'EngineMassRate': round(engine.G_air[0], 3),
        'CompTInSpecificHeat': round(engine.compressor.work_fluid._c_p_av_func(engine.compressor.T_in_stag), 3),
        'CompTOutSpecificHeat': round(engine.compressor.work_fluid._c_p_av_func(engine.compressor.output.T_out_stag[0]), 3),
        'SpHeatT0': 273,
        'CompTurbineTOutSpecificHeat':
            round(engine.compressor_turbine.work_fluid._c_p_av_func(engine.compressor_turbine.output.T_out_stag,
                                                                    alpha=engine.combustion_chamber.output.alpha)[0], 3),
        'CompTurbineTInSpecificHeat':
            round(engine.compressor_turbine.work_fluid._c_p_av_func(engine.compressor.T_in_stag,
                                                                    alpha=engine.combustion_chamber.output.alpha)[0], 3),
        'CompTurbineSpecificHeatLong': round(engine.compressor_turbine.work_fluid.c_p_av_int[0], 3),
        'FreeTurbineTOutSpecificHeat':
            round(engine.power_turbine.work_fluid._c_p_av_func(engine.power_turbine.output.T_out_stag,
                                                               alpha=engine.combustion_chamber.output.alpha)[0], 3),

        # Расчет на прочность
        'AngVelocity': strength_config.n,
        'BladeForce': round(strength_config.force_blade, 0),
        'TailWidth': round(FirstStageDisk().b_a_tail.value / 1e3, 5),
        'TailR1': round(FirstStageDisk().D_tail_in.value / 2 / 1e3, 5),
        'TailR2': round(FirstStageDisk().D_out.value / 2 / 1e3, 5),
        'OutTemp': strength_config.T_m,
        'BladeNumber': strength_config.z_blade,
        'Density': material.rho,
        'OutPressure': round(p_m / 1e6, 3),
        'FirstCalculationRows': get_first_calculation_rows(),
        'SecondCalculationRows': get_second_calculation_rows(),
        'RealStressRows': get_recalculated_stress_rows(),
        'AverageStressRows': get_average_stress_rows(),
        'SigmaRKFirst':  round(
            two_calculations.first_calc.sigma_r_arr[len(two_calculations.first_calc.sigma_r_arr) - 2][1] / 1e6, 3),
        'SigmaRKSecond': round(
            two_calculations.second_calc.sigma_r_arr[len(two_calculations.second_calc.sigma_r_arr) - 2][1] / 1e6, 3),
        'RecalculationCoef': round(two_calculations.k, 4),
        'MaxEqSigma': round(two_calculations.sigma_eq_max / 1e6, 3),
        'MinSafetyFactor': round(two_calculations.safety_factor_min, 3),
        'DHole1': round(strength_config.d_hole_arr[0], 4),
        'DHole2': round(strength_config.d_hole_arr[1], 4),
        'BHole1': round(strength_config.b_hole_arr[0], 4),
        'BHole2': round(strength_config.b_hole_arr[1], 4),
        'SigmaRHole1': round(two_calculations.sigma_r(FirstStageDisk().Dr_bolt.value / 1e3 / 2) / 1e6, 3),
        'SigmaTHole1': round(two_calculations.sigma_t(FirstStageDisk().Dr_bolt.value / 1e3 / 2) / 1e6, 3),
        'SigmaRHole2': round(two_calculations.sigma_r(0.0506) / 1e6, 3),
        'SigmaTHole2': round(two_calculations.sigma_t(0.0506) / 1e6, 3),
        'KHole1': round(two_calculations.k_hole_arr[0], 3),
        'KHole2': round(two_calculations.k_hole_arr[1], 3),
        'SigmaTKHole1': round(two_calculations.sigma_t_hole[0] / 1e6, 3),
        'SigmaTKHole2': round(two_calculations.sigma_t_hole[1] / 1e6, 3),

        # Расчет первой ступени
        'Reactivity': turbine[0].rho,
        'DeltaR': round(turbine.geom[0].delta_r_rk, 5),
        'StatorLRel': round(turbine.geom.l1_D1_ratio, 4),
        'StatorElongation': round(turbine.geom[0].l1_b_sa_ratio, 4),
        'RotorElongation': round(turbine.geom[0].l2_b_rk_ratio, 4),
        'StatorDeltaRel': round(turbine.geom[0].delta_a_b_sa_ratio, 4),
        'GammaIn': round(np.degrees(turbine.geom.gamma_in), 2),
        'GammaOut': round(np.degrees(turbine.geom.gamma_out), 2),
        'Ht': round(turbine[0].H0 / 1e6, 3),
        'Hc': round(turbine[0].H_s / 1e6, 3),
        'St1KGasShort': round(turbine[0].k_gas, 4),
        'Phi': round(turbine[0].phi, 3),
        'C1': round(turbine[0].c1, 3),
        'St1SpecificHeat': round(turbine[0].c_p_gas),
        'T1Prime': round(turbine[0].T1_ad),
        'St1PIn': round(turbine[0].p0_stag /1e6, 3),
        'P1': round(turbine[0].p1 / 1e6, 3),
        'St1GasTemp': round(turbine[0].T0_stag, 3),
        'St1GasConstant': round(turbine[0].work_fluid.R),
        'Rho1': round(turbine[0].rho1, 3),
        'Alpha1': round(np.degrees(turbine[0].alpha1), 3),
        'C1a': round(turbine[0].c1_a, 3),
        'G': round(turbine[0].G_turbine),
        'A1': round(turbine[0].A1_a, 5),
        'StatorDOut': round(turbine[0].D1, 4),
        'RotorDIn': round(turbine[0].D1, 4),
        'RotationSpeed': round(turbine[0].n, 3),
        'U1': round(turbine[0].u1, 3),
        'W1': round(turbine[0].w1, 3),
        'T1': round(turbine[0].T1, 3),
        'Hl': round(turbine[0].H_l / 1e6, 3),
        'BaRK': round(turbine.geom[0].b_rk, 3),
        'RotorDOut': round(turbine[0].D2, 4),
        'L2': round(turbine.geom[0].l2, 4),
        'RotorLRel': round(turbine.geom[0].l2 / turbine.geom[0].D2, 3),
        'U2': round(turbine[0].u2, 3),
        'Psi': round(turbine[0].psi, 3),
        'W2': round(turbine[0].w2, 3),
        'T2': round(turbine[0].T2, 3),
        'T2Prime': round(turbine[0].T2_ad, 3),
        'P2': round(turbine[0].p2 / 1e6, 3),
        'A2': round(turbine[0].A2_a, 6),
        'Rho2': round(turbine[0].rho2, 3),
        'C2a': round(turbine[0].c2_a, 3),
        'Beta2': round(np.degrees(turbine[0].beta2), 3),
        'Alpha2': round(np.degrees(turbine[0].alpha2), 3),
        'C2u': round(turbine[0].c2_u, 3),
        'C2': round(turbine[0].c2, 3),
        'C1u': round(turbine[0].c1_u, 3),
        'Lu': round(turbine[0].L_u / 1e6, 3),
        'EtaU': round(turbine[0].eta_u, 4),
        'hs': round(turbine[0].h_s / 1e3, 3),
        'hr': round(turbine[0].h_l / 1e3, 3),
        'hOut': round(turbine[0].h_v / 1e3, 3),
        'hRadial': round(turbine[0].h_z / 1e3, 3),
        'hVent': round(turbine[0].h_tv / 1e3, 6),
        'T2Stag': round(turbine[0].T_st_stag, 3),
        'P2Stag': round(turbine[0].p2_stag / 1e6, 3),
        'EtaPower': round(turbine[0].eta_t, 4),
        'Lt': round(turbine[0].L_t / 1e6, 3),
        'HtStag': round(turbine[0].H0_stag / 1e6, 3),
        'EtaT': round(turbine[0].eta_t_stag, 4),
        'St1KGasLong': round(turbine[0].work_fluid.k_av_int, 4),
        'St1KCalcError': round(turbine[0].dk_rel * 100, 3),
        'St1TInSpecificHeat': round(turbine[0].work_fluid._c_p_av_func(turbine[0].T0_stag,
                                                                       alpha=turbine[0].alpha_air), 3),
        'St1TOutSpecificHeat': round(turbine[0].work_fluid._c_p_av_func(turbine[0].T2, alpha=turbine[0].alpha_air), 3),
        'St1SpecificHeatLong': round(turbine[0].work_fluid.c_p_av_int, 3),


        # Расчет второй ступени
        'St2GasTemp': round(turbine[1].T0_stag, 3),
        'St2Reactivity': turbine[1].rho,
        'St2DeltaR': round(turbine.geom[1].delta_r_rk, 5),
        'St2L1': round(turbine.geom[1].l1, 4),
        'St2L2': round(turbine.geom[1].l2, 4),
        'St2Ht': round(turbine[1].H0 / 1e6, 3),
        'St2Hc': round(turbine[1].H_s / 1e6, 3),
        'St2KGasShort': round(turbine[1].k_gas, 4),
        'St2Phi': round(turbine[1].phi, 3),
        'St2C1': round(turbine[1].c1, 3),
        'St2SpecificHeat': round(turbine[1].c_p_gas),
        'St2T1Prime': round(turbine[1].T1_ad),
        'St2PIn': round(turbine[1].p0_stag / 1e6, 3),
        'St2P1': round(turbine[1].p1 / 1e6, 3),
        'St2GasConstant': round(turbine[1].work_fluid.R),
        'St2Rho1': round(turbine[1].rho1, 3),
        'St2Alpha1': round(np.degrees(turbine[1].alpha1), 3),
        'St2C1a': round(turbine[1].c1_a, 3),
        'St2G': round(turbine[1].G_turbine),
        'St2A1': round(turbine[1].A1_a, 5),
        'St2StatorDOut': round(turbine[1].D1, 4),
        'St2RotorDIn': round(turbine[1].D1, 4),
        'St2U1': round(turbine[1].u1, 3),
        'St2W1': round(turbine[1].w1, 3),
        'St2T1': round(turbine[1].T1, 3),
        'St2Hl': round(turbine[1].H_l / 1e6, 3),
        'St2BaRK': round(turbine.geom[1].b_rk, 3),
        'St2RotorDOut': round(turbine[1].D2, 4),
        'St2RotorLRel': round(turbine.geom[1].l2 / turbine.geom[0].D2, 3),
        'St2U2': round(turbine[1].u2, 3),
        'St2Psi': round(turbine[1].psi, 3),
        'St2W2': round(turbine[1].w2, 3),
        'St2T2': round(turbine[1].T2, 3),
        'St2T2Prime': round(turbine[1].T2_ad, 3),
        'St2P2': round(turbine[1].p2 / 1e6, 3),
        'St2A2': round(turbine[1].A2_a, 6),
        'St2Rho2': round(turbine[1].rho2, 3),
        'St2C2a': round(turbine[1].c2_a, 3),
        'St2Beta2': round(np.degrees(turbine[1].beta2), 3),
        'St2Alpha2': round(np.degrees(turbine[1].alpha2), 3),
        'St2C2u': round(turbine[1].c2_u, 3),
        'St2C2': round(turbine[1].c2, 3),
        'St2C1u': round(turbine[1].c1_u, 3),
        'St2Lu': round(turbine[1].L_u / 1e6, 3),
        'St2EtaU': round(turbine[1].eta_u, 4),
        'St2hs': round(turbine[1].h_s / 1e3, 3),
        'St2hr': round(turbine[1].h_l / 1e3, 3),
        'St2hOut': round(turbine[1].h_v / 1e3, 3),
        'St2hRadial': round(turbine[1].h_z / 1e3, 3),
        'St2hVent': round(turbine[1].h_tv / 1e3, 6),
        'St2T2Stag': round(turbine[1].T_st_stag, 3),
        'St2P2Stag': round(turbine[1].p2_stag / 1e6, 3),
        'St2EtaPower': round(turbine[1].eta_t, 4),
        'St2Lt': round(turbine[1].L_t / 1e6, 3),
        'St2HtStag': round(turbine[1].H0_stag / 1e6, 3),
        'St2EtaT': round(turbine[1].eta_t_stag, 4),
        'St2KGasLong': round(turbine[1].work_fluid.k_av_int, 4),
        'St2KCalcError': round(turbine[1].dk_rel * 100, 3),
        'St2TInSpecificHeat':
            round(turbine[1].work_fluid._c_p_av_func(turbine[1].T0_stag, alpha=turbine[1].alpha_air), 3),
        'St2TOutSpecificHeat':
            round(turbine[1].work_fluid._c_p_av_func(turbine[1].T2, alpha=turbine.alpha_air), 3),
        'St2SpecificHeatLong': round(turbine[1].work_fluid.c_p_av_int, 3),

        # Интегральные показатели турбины
        'TurbineLt': round(turbine.L_t_sum / 1e6, 3),
        'TurbineSpecificHeat': round(turbine.c_p_gas, 3),
        'TurbineKGas': round(turbine.k_gas, 3),
        'TurbineSpecificHeatStag': round(turbine.c_p_gas_stag, 3),
        'TurbineKGasStag': round(turbine.k_gas_stag, 3),
        'TurbineHt': round(turbine.H_t / 1e6, 3),
        'TurbineEtaT': round(turbine.eta_t, 4),
        'TurbineEtaL': round(turbine.eta_l, 4),
        'TurbineHtStag': round(turbine.H_t_stag / 1e6, 3),
        'TurbineEtaTStag': round(turbine.eta_t_stag, 4),
        'TurbinePower': round(turbine.N / 1e6, 3),
        'TurbineTOutStagSpecificHeat':
            round(turbine.work_fluid._c_p_av_func(turbine[1].T_st_stag, alpha=turbine.alpha_air), 3)
    }
    return result