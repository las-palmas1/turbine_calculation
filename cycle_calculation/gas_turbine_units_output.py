from gas import Air, KeroseneCombustionProducts
import numpy as np
import pandas


class CompressorOutput:
    def __init__(self):
        self._k_air = None
        self._p_out_stag = None
        self._T_out_stag = None
        self._eta_comp_stag = None
        self._c_p_air = None
        self._L_comp = None
        self._work_fluid = None
        self._output_frame = None

    @property
    def k_air(self) -> np.ndarray:
        return self._k_air

    @k_air.setter
    def k_air(self, value: np.ndarray):
        self._k_air = value

    @property
    def p_out_stag(self) -> np.ndarray:
        return self._p_out_stag

    @p_out_stag.setter
    def p_out_stag(self, value: np.ndarray):
        self._p_out_stag = value

    @property
    def T_out_stag(self) -> np.ndarray:
        return self._T_out_stag

    @T_out_stag.setter
    def T_out_stag(self, value: np.ndarray):
        self._T_out_stag = value

    @property
    def eta_comp_stag(self) -> np.ndarray:
        return self._eta_comp_stag

    @eta_comp_stag.setter
    def eta_comp_stag(self, value: np.ndarray):
        self._eta_comp_stag = value

    @property
    def c_p_air(self) -> np.ndarray:
        return self._c_p_air

    @c_p_air.setter
    def c_p_air(self, value: np.ndarray):
        self._c_p_air = value

    @property
    def L_comp(self) -> np.ndarray:
        return self._L_comp

    @L_comp.setter
    def L_comp(self, value: np.ndarray):
        self._L_comp = value

    @property
    def output_frame(self) -> pandas.DataFrame:
        return self._output_frame

    @output_frame.setter
    def output_frame(self, value: pandas.DataFrame):
        self._output_frame = value

    @property
    def work_fluid(self) -> Air:
        return self._work_fluid

    @work_fluid.setter
    def work_fluid(self, value: Air):
        self._work_fluid = value


class PowerTurbineOutput:
    def __init__(self):
        self._T_out_stag = None
        self._H_turb_stag = None
        self._eta_turb_stag = None
        self._k_gas = None
        self._c_p_gas = None
        self._L_turb = None
        self._pi_turb_stag = None
        self._output_frame = None
        self._work_fluid = None

    @property
    def pi_turb_stag(self) -> np.ndarray:
        return self._pi_turb_stag

    @pi_turb_stag.setter
    def pi_turb_stag(self, value: np.ndarray):
        self._pi_turb_stag = value

    @property
    def T_out_stag(self) -> np.ndarray:
        return self._T_out_stag

    @T_out_stag.setter
    def T_out_stag(self, value: np.ndarray):
        self._T_out_stag = value

    @property
    def H_turb_stag(self) -> np.ndarray:
        return self._H_turb_stag

    @H_turb_stag.setter
    def H_turb_stag(self, value: np.ndarray):
        self._H_turb_stag = value

    @property
    def eta_turb_stag(self) -> np.ndarray:
        return self._eta_turb_stag

    @eta_turb_stag.setter
    def eta_turb_stag(self, value: np.ndarray):
        self._eta_turb_stag = value

    @property
    def k_gas(self) -> np.ndarray:
        return self._k_gas

    @k_gas.setter
    def k_gas(self, value: np.ndarray):
        self._k_gas = value

    @property
    def c_p_gas(self) -> np.ndarray:
        return self._c_p_gas

    @c_p_gas.setter
    def c_p_gas(self, value: np.array):
        self._c_p_gas = value

    @property
    def L_turb(self) -> np.ndarray:
        return self._L_turb

    @L_turb.setter
    def L_turb(self, value: np.ndarray):
        self._L_turb = value

    @property
    def output_frame(self) -> pandas.DataFrame:
        return self._output_frame

    @output_frame.setter
    def output_frame(self, value: pandas.DataFrame):
        self._output_frame = value

    @property
    def work_fluid(self) -> KeroseneCombustionProducts:
        return self._work_fluid

    @work_fluid.setter
    def work_fluid(self, value: KeroseneCombustionProducts):
        self._work_fluid = value


class CompressorTurbineOutput:
    def __init__(self):
        self._L_turb = None
        self._T_out_stag = None
        self._c_p_gas = None
        self._k_gas = None
        self._pi_turb_stag = None
        self._p_out_stag = None
        self._H_turb_stag = None
        self._eta_turb_stag = None
        self._output_frame = None
        self._work_fluid = None

    @property
    def pi_turb_stag(self) -> np.ndarray:
        return self._pi_turb_stag

    @pi_turb_stag.setter
    def pi_turb_stag(self, value: np.ndarray):
        self._pi_turb_stag = value

    @property
    def T_out_stag(self) -> np.ndarray:
        return self._T_out_stag

    @T_out_stag.setter
    def T_out_stag(self, value: np.ndarray):
        self._T_out_stag = value

    @property
    def H_turb_stag(self) -> np.ndarray:
        return self._H_turb_stag

    @H_turb_stag.setter
    def H_turb_stag(self, value: np.ndarray):
        self._H_turb_stag = value

    @property
    def eta_turb_stag(self) -> np.ndarray:
        return self._eta_turb_stag

    @eta_turb_stag.setter
    def eta_turb_stag(self, value: np.ndarray):
        self._eta_turb_stag = value

    @property
    def k_gas(self) -> np.ndarray:
        return self._k_gas

    @k_gas.setter
    def k_gas(self, value: np.ndarray):
        self._k_gas = value

    @property
    def c_p_gas(self) -> np.ndarray:
        return self._c_p_gas

    @c_p_gas.setter
    def c_p_gas(self, value: np.array):
        self._c_p_gas = value

    @property
    def L_turb(self) -> np.ndarray:
        return self._L_turb

    @L_turb.setter
    def L_turb(self, value: np.ndarray):
        self._L_turb = value

    @property
    def p_out_stag(self) -> np.ndarray:
        return self._p_out_stag

    @p_out_stag.setter
    def p_out_stag(self, value: np.ndarray):
        self._p_out_stag = value

    @property
    def output_frame(self) -> pandas.DataFrame:
        return self._output_frame

    @output_frame.setter
    def output_frame(self, value: pandas.DataFrame):
        self._output_frame = value

    @property
    def work_fluid(self) -> KeroseneCombustionProducts:
        return self._work_fluid

    @work_fluid.setter
    def work_fluid(self, value: KeroseneCombustionProducts):
        self._work_fluid = value


class CombustionChamberOutput:
    def __init__(self):
        self._g_fuel = None
        self._g_gas = None
        self._alpha = None
        self._p_out_stag = None
        self._output_frame = None

    @property
    def g_fuel(self) -> np.ndarray:
        return self._g_fuel

    @g_fuel.setter
    def g_fuel(self, value: np.ndarray):
        self._g_fuel = value

    @property
    def g_gas(self) -> np.ndarray:
        return self._g_gas

    @g_gas.setter
    def g_gas(self, value: np.ndarray):
        self._g_gas = value

    @property
    def alpha(self) -> np.ndarray:
        return self._alpha

    @alpha.setter
    def alpha(self, value: np.ndarray):
        self._alpha = value

    @property
    def p_out_stag(self) -> np.ndarray:
        return self._p_out_stag

    @p_out_stag.setter
    def p_out_stag(self, value: np.ndarray):
        self._p_out_stag = value

    @property
    def output_frame(self) -> pandas.DataFrame:
        return self._output_frame

    @output_frame.setter
    def output_frame(self, value: pandas.DataFrame):
        self._output_frame = value


class InletOutput:
    def __init__(self):
        self._pi_v = None
        self._p_out_stag = None
        self._T_out_stag = None

    @property
    def pi_v(self):
        return self._pi_v

    @pi_v.setter
    def pi_v(self, value):
        self._pi_v = value

    @property
    def p_out_stag(self):
        return self._p_out_stag

    @p_out_stag.setter
    def p_out_stag(self, value):
        self._p_out_stag = value

    @property
    def T_out_stag(self):
        return self._T_out_stag

    @T_out_stag.setter
    def T_out_stag(self, value):
        self._T_out_stag = value


class OutletOutput:
    def __init__(self):
        self.T_out_stag = None
        self.p_in_stag = None
        self.T_out = None
        self.c_out = None
        self.p_out_stag = None
        self.k_gas = None
        self.c_p_gas = None
        self.output_frame = None
