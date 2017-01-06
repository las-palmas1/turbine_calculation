from gas import Air, KeroseneCombustionProducts
import numpy as np


class CompressorInput:
    def __init__(self):
        self._pi_comp_stag = None
        self._eta_comp_stag_p = None
        self._T_in_stag = None
        self._p_in_stag = None
        self._work_fluid = None

    @property
    def pi_comp_stag(self) -> np.ndarray:
        assert self._pi_comp_stag is not None, 'pi_comp_stag is None'
        return self._pi_comp_stag

    @pi_comp_stag.setter
    def pi_comp_stag(self, value: np.ndarray):
        self._pi_comp_stag = value

    @property
    def eta_comp_stag_p(self):
        assert self._eta_comp_stag_p is not None, 'eta_comp_stag_p is None'
        return self._eta_comp_stag_p

    @eta_comp_stag_p.setter
    def eta_comp_stag_p(self, value):
        self._eta_comp_stag_p = value

    @property
    def T_in_stag(self):
        assert self._T_in_stag is not None, 'T_in_stag is None'
        return self._T_in_stag

    @T_in_stag.setter
    def T_in_stag(self, value):
        self._T_in_stag = value

    @property
    def p_in_stag(self):
        assert self._p_in_stag is not None, 'p_in_stag is None'
        return self._p_in_stag

    @p_in_stag.setter
    def p_in_stag(self, value):
        self._p_in_stag = value

    @property
    def work_fluid(self) -> Air:
        assert self._work_fluid is not None, 'work_fluid is None'
        return self._work_fluid

    @work_fluid.setter
    def work_fluid(self, value: Air):
        self._work_fluid = value


class PowerTurbineInput:
    def __init__(self):
        self._p_in_stag = None
        self._p_out_stag = None
        self._T_in_stag = None
        self._eta_turb_stag_p = None
        self._alpha = None
        self._work_fluid = None
        self._T_out_stag_init = None

    @property
    def p_in_stag(self) -> np.ndarray:
        assert self._p_in_stag is not None, 'p_in_stag is None'
        return self._p_in_stag

    @p_in_stag.setter
    def p_in_stag(self, value: np.ndarray):
        self._p_in_stag = value

    @property
    def p_out_stag(self) -> np.ndarray:
        assert self._p_out_stag is not None, 'p_out is None'
        return self._p_out_stag

    @p_out_stag.setter
    def p_out_stag(self, value: np.ndarray):
        self._p_out_stag = value

    @property
    def T_in_stag(self) -> np.ndarray:
        assert self._T_in_stag is not None, 'T_in_stag is None'
        return self._T_in_stag

    @T_in_stag.setter
    def T_in_stag(self, value: np.ndarray):
        self._T_in_stag = value

    @property
    def eta_turb_stag_p(self):
        assert self._eta_turb_stag_p is not None, 'eta_turb_l is None'
        return self._eta_turb_stag_p

    @eta_turb_stag_p.setter
    def eta_turb_stag_p(self, value):
        self._eta_turb_stag_p = value

    @property
    def alpha(self) -> np.ndarray:
        assert self._alpha is not None, 'alpha is None'
        return self._alpha

    @alpha.setter
    def alpha(self, value: np.ndarray):
        self._alpha = value

    @property
    def work_fluid(self) -> KeroseneCombustionProducts:
        assert self._work_fluid is not None, 'work_fluid is None'
        return self._work_fluid

    @work_fluid.setter
    def work_fluid(self, value: KeroseneCombustionProducts):
        self._work_fluid = value

    @property
    def T_out_stag_init(self):
        assert self._T_out_stag_init is not None, 'T_out_stag_init is None'
        return self._T_out_stag_init

    @T_out_stag_init.setter
    def T_out_stag_init(self, value):
        self._T_out_stag_init = value


class CompressorTurbineInput:
    def __init__(self):
        self._g_gas = None
        self._L_comp = None
        self._T_in_stag = None
        self._p_in_stag = None
        self._eta_turb_stag_p = None
        self._alpha = None
        self._eta_m = None
        self._pi_turb_stag_init = None
        self._work_fluid = None

    @property
    def p_in_stag(self) -> np.ndarray:
        assert self._p_in_stag is not None, 'p_in_stag is None'
        return self._p_in_stag

    @p_in_stag.setter
    def p_in_stag(self, value: np.ndarray):
        self._p_in_stag = value

    @property
    def alpha(self) -> np.ndarray:
        assert self._alpha is not None, 'alpha is None'
        return self._alpha

    @alpha.setter
    def alpha(self, value: np.ndarray):
        self._alpha = value

    @property
    def g_gas(self) -> np.ndarray:
        assert self._g_gas is not None, 'g_gas is None'
        return self._g_gas

    @g_gas.setter
    def g_gas(self, value: np.ndarray):
        self._g_gas = value

    @property
    def L_comp(self) -> np.ndarray:
        assert self._L_comp is not None, 'L_comp is None'
        return self._L_comp

    @L_comp.setter
    def L_comp(self, value: np.ndarray):
        self._L_comp = value

    @property
    def T_in_stag(self) -> np.ndarray:
        assert self._T_in_stag is not None, 'T_in_stag is None'
        return self._T_in_stag

    @T_in_stag.setter
    def T_in_stag(self, value: np.ndarray):
        self._T_in_stag = value

    @property
    def eta_turb_stag_p(self):
        assert self._eta_turb_stag_p is not None, 'eta_turb_stag is None'
        return self._eta_turb_stag_p

    @eta_turb_stag_p.setter
    def eta_turb_stag_p(self, value):
        self._eta_turb_stag_p = value

    @property
    def eta_m(self):
        assert self._eta_m is not None, 'eta_m is None'
        return self._eta_m

    @eta_m.setter
    def eta_m(self, value):
        self._eta_m = value

    @property
    def work_fluid(self) -> KeroseneCombustionProducts:
        assert self._work_fluid is not None, 'work_fluid is None'
        return self._work_fluid

    @work_fluid.setter
    def work_fluid(self, value: KeroseneCombustionProducts):
        self._work_fluid = value

    @property
    def pi_turb_stag_init(self):
        assert self._pi_turb_stag_init is not None, 'pi_turb_stag_init is None'
        return self._pi_turb_stag_init

    @pi_turb_stag_init.setter
    def pi_turb_stag_init(self, value):
        self._pi_turb_stag_init = value


class CombustionChamberInput:
    def __init__(self):
        self._T_in_stag = None
        self._T_out_stag = None
        self._eta_comb = None
        self._Q_n = None
        self._l0 = None
        self._p_in_stag = None
        self._sigma_comb = None
        self._work_fluid_in = None
        self._work_fluid_out = None
        self._g_outflow = None
        self._g_cooling = None
        self._g_return = None

    @property
    def T_in_stag(self) -> np.ndarray:
        assert self._T_in_stag is not None, 'T_in_stag is None'
        return self._T_in_stag

    @T_in_stag.setter
    def T_in_stag(self, value: np.ndarray):
        self._T_in_stag = value

    @property
    def p_in_stag(self) -> np.ndarray:
        assert self._p_in_stag is not None, 'p_in_stag is None'
        return self._p_in_stag

    @p_in_stag.setter
    def p_in_stag(self, value: np.ndarray):
        self._p_in_stag = value

    @property
    def T_out_stag(self):
        assert self._T_out_stag is not None, 'T_out_stag is None'
        return self._T_out_stag

    @T_out_stag.setter
    def T_out_stag(self, value):
        self._T_out_stag = value

    @property
    def eta_comb(self):
        assert self._eta_comb is not None, 'eta_g is None'
        return self._eta_comb

    @eta_comb.setter
    def eta_comb(self, value):
        self._eta_comb = value

    @property
    def sigma_comb(self):
        assert self._sigma_comb is not None, 'sigma_comb is None'
        return self._sigma_comb

    @sigma_comb.setter
    def sigma_comb(self, value):
        self._sigma_comb = value

    @property
    def Q_n(self):
        assert self._Q_n is not None, 'Q_n is None'
        return self._Q_n

    @Q_n.setter
    def Q_n(self, value):
        self._Q_n = value

    @property
    def l0(self):
        assert self._l0 is not None, 'l0 is None'
        return self._l0

    @l0.setter
    def l0(self, value):
        self._l0 = value

    @property
    def work_fluid_in(self) -> Air:
        assert self._work_fluid_in is not None, 'work_fluid_in is None'
        return self._work_fluid_in

    @work_fluid_in.setter
    def work_fluid_in(self, value: Air):
        self._work_fluid_in = value

    @property
    def work_fluid_out(self) -> KeroseneCombustionProducts:
        assert self._work_fluid_out is not None, 'work_fluid_out is None'
        return self._work_fluid_out

    @work_fluid_out.setter
    def work_fluid_out(self, value: KeroseneCombustionProducts):
        self._work_fluid_out = value

    @property
    def g_outflow(self):
        assert self._g_outflow is not None, 'g_outflow is None'
        return self._g_outflow

    @g_outflow.setter
    def g_outflow(self, value):
        self._g_outflow = value

    @property
    def g_cooling(self):
        assert self._g_cooling is not None, 'g_cooling is None'
        return self._g_cooling

    @g_cooling.setter
    def g_cooling(self, value):
        self._g_cooling = value

    @property
    def g_return(self):
        assert self._g_return is not None, 'g_return is None'
        return self._g_return

    @g_return.setter
    def g_return(self, value):
        self._g_return = value


class InletInput:
    def __init__(self):
        self._T_in = None
        self._p_in = None
        self._M_v = None
        self._v = None
        self._sigma_in = None
        self._work_fluid = None

    @property
    def T_in(self):
        assert self._T_in is not None, 'T_in is None'
        return self._T_in

    @T_in.setter
    def T_in(self, value):
        self._T_in = value

    @property
    def p_in(self):
        assert self._p_in is not None, 'p_in is None'
        return self._p_in

    @p_in.setter
    def p_in(self, value):
        self._p_in = value

    @property
    def M_v(self):
        return self._M_v

    @M_v.setter
    def M_v(self, value):
        self._M_v = value

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, value):
        self._v = value

    @property
    def sigma_in(self):
        assert self._sigma_in is not None, 'sigma_in is None'
        return self._sigma_in

    @sigma_in.setter
    def sigma_in(self, value):
        self._sigma_in = value

    @property
    def work_fluid(self) -> Air:
        assert self._work_fluid is not None, 'work_fluid is None'
        return self._work_fluid

    @work_fluid.setter
    def work_fluid(self, value: Air):
        self._work_fluid = value


class OutletInput:
    def __init__(self):
        self._work_fluid = None
        self._alpha = None
        self._p_out = None
        self._T_in_stag = None
        self._lam_out = None
        self._sigma_out = None

    @property
    def work_fluid(self) -> KeroseneCombustionProducts:
        assert self._work_fluid is not None, 'work_fluid is None'
        return self._work_fluid

    @work_fluid.setter
    def work_fluid(self, value: KeroseneCombustionProducts):
        self._work_fluid = value

    @property
    def alpha(self):
        assert self._alpha is not None, 'alpha is None'
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value

    @property
    def p_out(self):
        assert self._p_out is not None, 'p_out is None'
        return self._p_out

    @p_out.setter
    def p_out(self, value):
        self._p_out = value

    @property
    def lam_out(self):
        assert self._lam_out is not None, 'lam_out is None'
        return self._lam_out

    @lam_out.setter
    def lam_out(self, value):
        self._lam_out = value

    @property
    def sigma_out(self):
        assert self._sigma_out is not None, 'sigma_out is None'
        return self._sigma_out

    @sigma_out.setter
    def sigma_out(self, value):
        self._sigma_out = value

    @property
    def T_in_stag(self) -> np.ndarray:
        assert self._T_in_stag is not None, 'T_in_stag is None'
        return self._T_in_stag

    @T_in_stag.setter
    def T_in_stag(self, value: np.ndarray):
        self._T_in_stag = value