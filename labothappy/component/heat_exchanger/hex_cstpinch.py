import __init__

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

from CoolProp.CoolProp import AbstractState
import CoolProp.CoolProp as CoolProp
from scipy.optimize import fsolve, root, minimize, brentq
import numpy as np
import math

# Fluids not available in CoolProp — properties must come from the connector directly
CUSTOM_FLUIDS = {'SolarSalt'}


def _get_abstract_state(fluid):
    """Return CoolProp AbstractState, or None for custom fluids."""
    if fluid in CUSTOM_FLUIDS:
        return None
    return AbstractState("HEOS", fluid)


def _h_from_p_T(AS, fluid, p, T, su_connector=None):
    """Enthalpy [J/kg] from (p, T). Uses custom correlations for non-CoolProp fluids."""
    if fluid in CUSTOM_FLUIDS:
        return type(su_connector)._h(T)
    AS.update(CoolProp.PT_INPUTS, p, T)
    return AS.hmass()


def _T_from_h_p(AS, fluid, h, p, su_connector=None):
    """Temperature [K] from (h, p). Uses custom correlations for non-CoolProp fluids."""
    if fluid in CUSTOM_FLUIDS:
        return type(su_connector)._T_from_h(h)
    AS.update(CoolProp.HmassP_INPUTS, h, p)
    return AS.T()


def _cp_from_p_T(AS, fluid, p, T, su_connector=None):
    """Specific heat [J/(kg·K)] from (p, T). Uses custom correlations for non-CoolProp fluids."""
    if fluid in CUSTOM_FLUIDS:
        return type(su_connector)._Cp(T)
    AS.update(CoolProp.PT_INPUTS, p, T)
    return AS.cpmass()


class HexCstPinch(BaseComponent):
    """
    Component: Heat Exchanger with constant pinch point.

    Adapted to support custom fluids (e.g. SolarSalt) that bypass CoolProp.
    When the hot-side fluid is a custom fluid, thermophysical properties are
    computed from the connector's own static correlations instead of CoolProp.

    Note: phase-change logic (evaporator / condenser) is only supported for
    standard CoolProp fluids. Custom fluids are treated as single-phase.
    """

    def __init__(self):
        super().__init__()
        self.su_C = MassConnector()
        self.su_H = MassConnector()
        self.ex_C = MassConnector()
        self.ex_H = MassConnector()

        self.DP_c = 0
        self.DP_h = 0

        self.Q = HeatConnector()
        self.guesses = {}

    def get_required_inputs(self):
        return ['fluid_C', 'h_su_C', 'P_su_C', 'm_dot_C', 'fluid_H', 'h_su_H', 'P_su_H', 'm_dot_H']

    def get_required_parameters(self):
        return ['Pinch', 'Delta_T_sh_sc', 'HX_type']

    def get_required_guesses(self):
        return ['P_sat']

    # ------------------------------------------------------------------
    # Single-phase solver for custom fluids (e.g. SolarSalt hot side)
    # ------------------------------------------------------------------

    def _solve_single_phase_custom_H(self):
        """
        Simplified single-phase solver used when the hot-side fluid is custom.
        Uses Cp-based capacity rates and enforces the pinch point at the outlet.

        For a counter-flow heat exchanger with a custom hot fluid (no phase
        change), the pinch always occurs at one of the two ends.  We set Q_dot
        such that min(T_H_outlet - T_C_inlet, T_H_inlet - T_C_outlet) = Pinch.
        """
        cp_H = type(self.su_H)._Cp(self.su_H.T)
        C_H  = self.su_H.m_dot * cp_H

        # Cold side Cp (CoolProp)
        cp_C = _cp_from_p_T(self.AS_C, self.su_C.fluid,
                             self.su_C.p, self.su_C.T,
                             su_connector=self.su_C)
        C_C  = self.su_C.m_dot * cp_C

        C_min = min(C_H, C_C)
        C_max = max(C_H, C_C)

        # Q_dot constrained by pinch at the hot outlet (most common for C_H > C_C)
        # T_H_out = T_H_in - Q / C_H  >=  T_C_in + Pinch
        # Q_pinch_hot_outlet = C_H * (T_H_in - T_C_in - Pinch)

        Q_pinch_H_out = C_H * (self.su_H.T - self.su_C.T - self.params['Pinch'])
        Q_pinch_C_out = C_C * (self.su_H.T - self.params['Pinch'] - self.su_C.T)

        Q_dot = min(Q_pinch_H_out, Q_pinch_C_out,
                    C_min * (self.su_H.T - self.su_C.T))
        Q_dot = max(Q_dot, 0.0)

        self.Q_dot_val = Q_dot

        # Outlet enthalpies
        h_ex_H = self.su_H.h - Q_dot / self.su_H.m_dot
        h_ex_C = self.su_C.h + Q_dot / self.su_C.m_dot

        # Outlet temperatures
        T_ex_H = type(self.su_H)._T_from_h(h_ex_H)

        self.AS_C.update(CoolProp.HmassP_INPUTS, h_ex_C, self.su_C.p)
        T_ex_C = self.AS_C.T()

        # Update connectors — hot side (custom fluid, no CoolProp)
        self.ex_H.fluid  = self.su_H.fluid
        self.ex_H.m_dot  = self.su_H.m_dot
        self.ex_H.p      = self.su_H.p - self.DP_h
        self.ex_H.h      = h_ex_H
        self.ex_H.T      = T_ex_H
        self.ex_H.cp     = type(self.su_H)._Cp(T_ex_H)
        self.ex_H.D      = type(self.su_H)._rho(T_ex_H)
        self.ex_H.state_known      = True
        self.ex_H.completely_known = True

        # Cold side (standard CoolProp fluid)
        self.ex_C.set_fluid(self.su_C.fluid)
        self.ex_C.set_m_dot(self.su_C.m_dot)
        self.ex_C.set_h(h_ex_C)
        self.ex_C.set_p(self.su_C.p - self.DP_c)

        self.Q.set_Q_dot(Q_dot)
        self.solved = True

    # ------------------------------------------------------------------
    # Original pinch-based solver (standard CoolProp fluids only)
    # ------------------------------------------------------------------

    def equivalent_effectiveness(self):
        if self.params['HX_type'] == "evaporator":
            self.AS_H.update(CoolProp.PT_INPUTS, self.su_H.p, self.su_C.T)
            H_h_id = self.AS_H.hmass()
            self.AS_C.update(CoolProp.PT_INPUTS, self.P_sat, self.su_H.T)
            H_c_id = self.AS_C.hmass()
        elif self.params['HX_type'] == "condenser":
            self.AS_H.update(CoolProp.PT_INPUTS, self.P_sat, self.su_C.T)
            H_h_id = self.AS_H.hmass()
            self.AS_C.update(CoolProp.PT_INPUTS, self.su_C.p, self.su_H.T)
            H_c_id = self.AS_C.hmass()

        Q_dot_maxh = self.su_H.m_dot * (self.su_H.h - H_h_id)
        Q_dot_maxc = self.su_C.m_dot * (H_c_id - self.su_C.h)
        self.Q_dot_max_ext = np.min([Q_dot_maxh, Q_dot_maxc])

        Q_dot_save  = self.Q_dot
        Pinch_save  = self.params['Pinch']
        P_sat_save  = self.P_sat
        su_C_p_save = self.su_C.p; su_H_p_save = self.su_H.p
        su_C_T_save = self.su_C.T; su_H_T_save = self.su_H.T
        su_C_h_save = self.su_C.h; su_H_h_save = self.su_H.h

        self.params['Pinch'] = 1e-2

        if self.params['HX_type'] == 'evaporator':
            guess_T_sat = self.su_H.T - self.params['Pinch'] - self.params['Delta_T_sh_sc']
            self.AS_C.update(CoolProp.QT_INPUTS, 0.5, guess_T_sat)
            P_ev_guess = self.AS_C.p()
            try:
                root(self.system_evap, [P_ev_guess], method='lm', tol=1e-7)
            except Exception as e:
                print(f"Convergence problem in pinch analysis of evaporator model: {e}")
        elif self.params['HX_type'] == 'condenser':
            guess_T_sat = self.su_C.T + self.params['Pinch'] + self.params['Delta_T_sh_sc']
            self.AS_H.update(CoolProp.QT_INPUTS, 0.5, guess_T_sat)
            P_cd_guess = self.AS_H.p()
            try:
                fsolve(self.system_cond, [P_cd_guess])
            except Exception as e:
                print(f"Convergence problem in pinch analysis of condenser model: {e}")

        self.Q_dot_max_int = self.Q_dot
        self.params['Pinch'] = Pinch_save
        self.su_C.set_p(su_C_p_save); self.su_H.set_p(su_H_p_save)
        self.su_C.set_T(su_C_T_save); self.su_H.set_T(su_H_T_save)
        self.su_C.h = su_C_h_save;    self.su_H.h = su_H_h_save
        self.P_sat  = P_sat_save

        self.Q_dot_max = np.min([self.Q_dot_max_ext, self.Q_dot_max_int])
        self.Q_dot = Q_dot_save

        if np.isfinite(self.Q_dot_max) and self.Q_dot_max > 0:
            self.epsilon = self.Q_dot / self.Q_dot_max
        else:
            self.epsilon = np.nan

    def system_evap(self, x):
        self.Q_dot_sc = 0; self.Q_dot_tp = 0; self.Q_dot_sh = 0
        self.P_ev = x
        PP_list = []
        self.PP_array = []

        P_triple = self.AS_C.trivial_keyed_output(CoolProp.iP_triple)
        if self.P_ev < P_triple:
            self.P_ev = P_triple + 1 - self.DP_c
        P_crit = self.AS_C.trivial_keyed_output(CoolProp.iP_critical)
        if self.P_ev - self.DP_c > P_crit:
            self.P_ev = P_crit - 1000 - self.DP_c

        self.x_start_evap = max(0, self.su_C.x)

        self.P_ev_x0 = self.P_ev + self.DP_c / 2
        self.P_ev_x1 = self.P_ev - self.DP_c / 2

        self.AS_C.update(CoolProp.PQ_INPUTS, self.P_ev_x0, self.x_start_evap)
        T_ev_x0 = self.AS_C.T(); self.T_ev_x0 = T_ev_x0
        self.AS_C.update(CoolProp.PQ_INPUTS, self.P_ev_x1, 1)
        T_ev_x1 = self.AS_C.T(); self.T_ev_x1 = T_ev_x1

        self.h_C_su = self.su_C.h
        self.AS_C.update(CoolProp.PQ_INPUTS, self.P_ev_x0, self.x_start_evap)
        self.h_C_x0 = self.AS_C.hmass()
        self.Q_dot_sc = self.su_C.m_dot * max(self.h_C_x0 - self.h_C_su, 0)

        self.AS_C.update(CoolProp.PQ_INPUTS, self.P_ev_x1, 1)
        h_C_x1 = self.AS_C.hmass(); self.h_C_x1 = h_C_x1
        self.Q_dot_tp = self.su_C.m_dot * min(self.h_C_x1 - self.h_C_x0,
                                               self.h_C_x1 - self.su_C.h)

        self.T_C_ex = self.T_ev_x1 + self.params['Delta_T_sh_sc']
        self.AS_C.update(CoolProp.PT_INPUTS, self.P_ev_x1, self.T_C_ex)
        if self.AS_C.cpmass() < 0:
            h_C_ex = CoolProp.PropsSI('H', 'P', self.P_ev_x1, 'T', self.T_C_ex, self.su_C.fluid)
        else:
            h_C_ex = self.AS_C.hmass()
            if h_C_ex < h_C_x1:
                h_C_ex = CoolProp.PropsSI('H', 'P', self.P_ev_x1, 'T', self.T_C_ex, self.su_C.fluid)
        self.h_C_ex = h_C_ex
        delta_h_sh = self.h_C_ex - max(self.h_C_x1, self.su_C.h)
        self.Q_dot_sh = self.su_C.m_dot * max(0.0, delta_h_sh)

        Q_dot_ev = self.Q_dot_sc + self.Q_dot_tp + self.Q_dot_sh

        self.h_H_x1 = self.su_H.h - self.Q_dot_sh / self.su_H.m_dot
        self.AS_H.update(CoolProp.HmassP_INPUTS, self.h_H_x1, self.su_H.p)
        self.T_H_x1 = self.AS_H.T()

        self.h_H_x0 = self.h_H_x1 - self.Q_dot_tp / self.su_H.m_dot
        self.AS_H.update(CoolProp.HmassP_INPUTS, self.h_H_x0, self.su_H.p)
        self.T_H_x0 = self.AS_H.T()

        self.h_H_ex = self.h_H_x0 - self.Q_dot_sc / self.su_H.m_dot
        self.AS_H.update(CoolProp.HmassP_INPUTS, self.h_H_ex, self.su_H.p)
        self.T_H_ex = self.AS_H.T()

        PP_list.append(self.T_H_ex - self.su_C.T)
        if self.Q_dot_sc > 0:
            PP_list.append(self.T_H_x0 - self.T_ev_x0)
        if self.Q_dot_sh > 0:
            PP_list.append(self.su_H.T - self.T_C_ex)

        self.PP_array = np.array(PP_list)
        self.PPTD = min(self.PP_array)
        self.res = self.PPTD - self.params['Pinch']
        self.Q_dot = Q_dot_ev
        self.P_sat = self.P_ev
        return self.res

    def system_cond(self, x):
        self.Q_dot_sc = 0; self.Q_dot_tp = 0; self.Q_dot_sh = 0
        P_cd = x; self.P_cd = P_cd
        PP_list = []; self.PP_array = []

        P_triple = self.AS_H.trivial_keyed_output(CoolProp.iP_triple)
        if self.P_cd < P_triple:
            self.P_cd = P_triple + 1 - self.DP_h
        P_crit = self.AS_H.trivial_keyed_output(CoolProp.iP_critical)
        if self.P_cd - self.DP_h > P_crit:
            self.P_cd = P_crit - 1000 - self.DP_h

        self.x_start_cd = min(1, abs(self.su_H.x))

        self.P_cd_x1 = self.P_cd + self.DP_h / 2
        self.P_cd_x0 = self.P_cd - self.DP_h / 2

        self.AS_H.update(CoolProp.PQ_INPUTS, self.P_cd_x1, self.x_start_cd)
        T_cd_x1 = self.AS_H.T(); self.h_H_x1 = self.AS_H.hmass(); self.T_cd_x1 = T_cd_x1
        self.AS_H.update(CoolProp.PQ_INPUTS, self.P_cd_x0, 0)
        T_cd_x0 = self.AS_H.T(); h_H_x0 = self.AS_H.hmass(); self.T_cd_x0 = T_cd_x0

        self.h_H_su = self.su_H.h
        self.Q_dot_sh = self.su_H.m_dot * max(self.h_H_su - self.h_H_x1, 0)
        self.h_H_x0 = h_H_x0
        self.Q_dot_tp = self.su_H.m_dot * max(
            min(self.h_H_x1 - self.h_H_x0, self.su_H.h - self.h_H_x0), 0)

        self.T_H_ex = self.T_cd_x0 - self.params['Delta_T_sh_sc']
        self.AS_H.update(CoolProp.PT_INPUTS, self.P_cd_x0, self.T_H_ex)
        if self.AS_H.cpmass() < 0:
            h_H_ex = CoolProp.PropsSI('H', 'P', self.P_cd_x0, 'T', self.T_H_ex, self.su_H.fluid)
        else:
            h_H_ex = self.AS_H.hmass()
            if h_H_ex > h_H_x0:
                h_H_ex = CoolProp.PropsSI('H', 'P', self.P_cd_x0, 'T', self.T_H_ex, self.su_H.fluid)
        self.h_H_ex = h_H_ex
        delta_h_sc = min(self.h_H_x0, self.su_H.h) - self.h_H_ex
        self.Q_dot_sc = self.su_H.m_dot * max(0.0, delta_h_sc)

        Q_dot_cd = self.Q_dot_sc + self.Q_dot_tp + self.Q_dot_sh

        self.h_C_x0 = self.su_C.h + self.Q_dot_sc / self.su_C.m_dot
        self.AS_C.update(CoolProp.HmassP_INPUTS, self.h_C_x0, self.su_C.p)
        self.T_C_x0 = self.AS_C.T()

        self.h_C_x1 = self.h_C_x0 + self.Q_dot_tp / self.su_C.m_dot
        self.AS_C.update(CoolProp.HmassP_INPUTS, self.h_C_x1, self.su_C.p)
        self.T_C_x1 = self.AS_C.T()

        self.h_C_ex = self.h_C_x1 + self.Q_dot_sh / self.su_C.m_dot
        self.AS_C.update(CoolProp.HmassP_INPUTS, self.h_C_ex, self.su_C.p)
        self.T_C_ex = self.AS_C.T()

        PP_list.append(self.su_H.T - self.T_C_ex)
        if self.Q_dot_tp > 0:
            PP_list.append(self.T_cd_x1 - self.T_C_x1)
        if self.Q_dot_sc > 0:
            PP_list.append(self.T_H_ex - self.su_C.T)

        self.PP_array = np.array(PP_list)
        self.PPTD = np.min(self.PP_array)
        self.res   = self.PPTD - self.params['Pinch']
        self.Q_dot = Q_dot_cd
        self.P_sat = P_cd
        return self.res

    def solve(self):
        self.change_flag = 0
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            print("HTX IS NOT CALCULABLE")
            return

        self.DP_h = self.params.get('DP_h', 0)
        self.DP_c = self.params.get('DP_c', 0)

        # Build AbstractState — None for custom fluids
        self.AS_C = _get_abstract_state(self.su_C.fluid)
        self.AS_H = _get_abstract_state(self.su_H.fluid)

        # --- Custom hot-side fluid path (e.g. SolarSalt) ---
        if self.su_H.fluid in CUSTOM_FLUIDS:
            self._solve_single_phase_custom_H()
            return

        # --- Standard CoolProp path (unchanged from original) ---
        if self.params['HX_type'] == 'evaporator':
            guess_T_sat_max = self.su_H.T
            self.AS_C.update(CoolProp.QT_INPUTS, 0.5, guess_T_sat_max)
            P_ev_guess = self.AS_C.p()

            P_triple = self.AS_C.trivial_keyed_output(CoolProp.iP_triple)
            P_crit   = self.AS_C.trivial_keyed_output(CoolProp.iP_critical)

            max_iter = 1000; step = 0.1; max_step = 5.0
            lower_bound = max(P_triple * 1.1 + self.DP_c / 2, 0.5 * P_ev_guess)
            upper_bound = P_ev_guess

            for _ in range(max_iter):
                res1 = self.system_evap(lower_bound)
                res2 = self.system_evap(upper_bound)
                if res1 * res2 <= 0:
                    break
                self.su_C.set_T(self.su_C.T - step)
                step = min(step * 2, max_step)

            if res1 * res2 > 0:
                raise ValueError("Could not bracket evaporator pressure root.")

            self.P_solution, self.results = brentq(
                self.system_evap,
                max(P_triple * 1.1 + self.DP_c / 2, 0.1 * P_ev_guess),
                upper_bound, xtol=1e-6, rtol=1e-8, maxiter=100, full_output=True)
            self.system_evap(self.P_solution)
            self.update_connectors()
            self.solved = self.results.converged and abs(self.res) < 1e-1

        elif self.params['HX_type'] == 'condenser':
            guess_T_sat_min = self.su_C.T
            self.AS_H.update(CoolProp.QT_INPUTS, 0.5, guess_T_sat_min)
            P_cd_guess = self.AS_H.p()

            P_triple = self.AS_H.trivial_keyed_output(CoolProp.iP_triple)
            P_crit   = self.AS_H.trivial_keyed_output(CoolProp.iP_critical)

            max_iter = 100; step = 0.1; max_step = 5.0
            for _ in range(max_iter):
                self.AS_H.update(CoolProp.QT_INPUTS, 0.5, self.su_C.T)
                lower_bound = self.AS_H.p()
                self.AS_H.update(CoolProp.QT_INPUTS, 0.5, self.su_H.T)
                upper_bound = max(P_crit * 0.95, self.AS_H.p())
                res1 = self.system_cond(lower_bound)
                res2 = self.system_cond(upper_bound)
                if res1 * res2 <= 0:
                    break
                self.su_H.set_T(self.su_H.T + step)
                step = min(step * 2, max_step)

            if res1 * res2 > 0:
                raise ValueError("Could not bracket condenser pressure root.")

            self.P_solution, self.results = brentq(
                self.system_cond,
                max(P_triple * 1.1 + self.DP_h, 0.1 * P_cd_guess),
                upper_bound, xtol=1e-6, rtol=1e-8, maxiter=100, full_output=True)
            self.system_cond(self.P_solution)
            self.update_connectors()
            self.solved = self.results.converged and abs(self.res) < 1e-1

    def update_connectors(self):
        if self.params['HX_type'] == 'evaporator':
            self.su_C.set_p(self.P_ev_x0)
            self.ex_C.set_fluid(self.su_C.fluid)
            self.ex_C.set_p(self.P_ev_x1)
            self.ex_C.set_h(self.h_C_ex)
            self.ex_C.set_m_dot(self.su_C.m_dot)

            # Hot side
            if self.su_H.fluid in CUSTOM_FLUIDS:
                h_ex_H = self.h_H_ex
                self.ex_H.fluid = self.su_H.fluid
                self.ex_H.m_dot = self.su_H.m_dot
                self.ex_H.p     = self.su_H.p - self.DP_h
                self.ex_H.h     = h_ex_H
                self.ex_H.T     = type(self.su_H)._T_from_h(h_ex_H)
                self.ex_H.cp    = type(self.su_H)._Cp(self.ex_H.T)
                self.ex_H.D     = type(self.su_H)._rho(self.ex_H.T)
                self.ex_H.state_known = True; self.ex_H.completely_known = True
            else:
                self.ex_H.set_fluid(self.su_H.fluid)
                self.ex_H.set_m_dot(self.su_H.m_dot)
                self.ex_H.set_p(self.su_H.p - self.DP_h)
                self.ex_H.set_T(self.T_H_ex)

            self.Q.set_Q_dot(self.Q_dot)

        else:  # condenser
            self.su_H.set_p(self.P_cd_x1)
            self.ex_H.set_fluid(self.su_H.fluid)
            self.ex_H.set_p(self.P_cd_x0)
            self.ex_H.set_h(self.h_H_ex)
            self.ex_H.set_m_dot(self.su_H.m_dot)

            self.ex_C.set_fluid(self.su_C.fluid)
            self.ex_C.set_m_dot(self.su_C.m_dot)
            self.ex_C.set_p(self.su_C.p - self.DP_c)
            self.ex_C.set_T(self.T_C_ex)

            self.Q.set_Q_dot(self.Q_dot)

    def print_results(self):
        print("=== Heat Exchanger Results ===")
        print(f"Q: {self.Q.Q_dot}")
        if self.params['HX_type'] == 'evaporator':
            print(f"P_sat: {self.su_C.p}")
        else:
            print(f"P_sat: {self.su_H.p}")
        print("======================")

    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, P={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, P={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, P={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, P={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q.Q_dot}")
        print("======================")

    def plot_disc(self):
        import matplotlib.pyplot as plt
        if self.params['HX_type'] == 'condenser':
            plt.figure()
            plt.plot([0, self.Q_dot_sh], [self.su_H.T, self.T_cd_x1], 'r', label='H')
            plt.plot([self.Q_dot_sh, self.Q_dot_sh + self.Q_dot_tp], [self.T_cd_x1, self.T_cd_x0], 'r')
            plt.plot([self.Q_dot_sh + self.Q_dot_tp, self.Q_dot], [self.T_cd_x0, self.ex_H.T], 'r')
            plt.plot([0, self.Q_dot_sh], [self.ex_C.T, self.T_C_x1], 'b', label='C')
            plt.plot([self.Q_dot_sh, self.Q_dot_sh + self.Q_dot_tp], [self.T_C_x1, self.T_C_x0], 'b')
            plt.plot([self.Q_dot_sh + self.Q_dot_tp, self.Q_dot], [self.T_C_x0, self.su_C.T], 'b')
            plt.grid(); plt.legend(); plt.show()
        if self.params['HX_type'] == 'evaporator':
            plt.figure()
            plt.plot([0, self.Q_dot_sc], [self.su_C.T, self.T_ev_x0], 'b', label='C')
            plt.plot([self.Q_dot_sc, self.Q_dot_sc + self.Q_dot_tp], [self.T_ev_x0, self.T_ev_x1], 'b')
            plt.plot([self.Q_dot_sc + self.Q_dot_tp, self.Q_dot], [self.T_ev_x1, self.ex_C.T], 'b')
            plt.plot([0, self.Q_dot_sc], [self.ex_H.T, self.T_H_x0], 'r', label='H')
            plt.plot([self.Q_dot_sc, self.Q_dot_sc + self.Q_dot_tp], [self.T_H_x0, self.T_H_x1], 'r')
            plt.plot([self.Q_dot_sc + self.Q_dot_tp, self.Q_dot], [self.T_H_x1, self.su_H.T], 'r')
            plt.grid(); plt.legend(); plt.show()