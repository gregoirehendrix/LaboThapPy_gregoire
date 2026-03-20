from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

import CoolProp.CoolProp as CP
import numpy as np
import math

# Fluids not available in CoolProp — properties must come from the connector directly
CUSTOM_FLUIDS = {'SolarSalt'}


def _get_abstract_state(fluid):
    """
    Return a CoolProp AbstractState for standard fluids.
    Returns None for custom fluids (e.g. SolarSalt) that bypass CoolProp.
    """
    if fluid in CUSTOM_FLUIDS:
        return None
    return CP.AbstractState('BICUBIC&HEOS', fluid)


def _get_T_from_h_p(AS, fluid, h, p, su_connector=None):
    """
    Return temperature [K] from enthalpy and pressure.
    For custom fluids, uses the connector's own _T_from_h method.
    """
    if fluid in CUSTOM_FLUIDS:
        # Recover T from h using the custom connector's static correlation
        connector_class = type(su_connector)
        return connector_class._T_from_h(h)
    AS.update(CP.HmassP_INPUTS, h, p)
    return AS.T()


def _get_h_from_p_T(AS, fluid, p, T, su_connector=None):
    """
    Return enthalpy [J/kg] from pressure and temperature.
    For custom fluids, uses the connector's own _h method.
    """
    if fluid in CUSTOM_FLUIDS:
        connector_class = type(su_connector)
        return connector_class._h(T)
    AS.update(CP.PT_INPUTS, p, T)
    return AS.hmass()


def _get_cp_from_p_T(AS, fluid, p, T, su_connector=None):
    """
    Return specific heat capacity [J/(kg·K)] from pressure and temperature.
    For custom fluids, uses the connector's own _Cp method.
    """
    if fluid in CUSTOM_FLUIDS:
        connector_class = type(su_connector)
        return connector_class._Cp(T)
    AS.update(CP.PT_INPUTS, p, T)
    return AS.cpmass()


class HexCstEffDisc(BaseComponent):
    """
    Component: Counterflow Heat Exchanger with Constant Effectiveness (HXEffCstDisc)
    
    Model: Discretized Counterflow Heat Exchanger with Fixed Effectiveness and Pinch Check
    
    Adapted to support custom fluids (e.g. SolarSalt) that bypass CoolProp.
    For custom fluids, thermophysical properties are computed from the connector's
    own static correlations instead of CoolProp AbstractState.
    """
    
    def __init__(self):
        super().__init__()
        self.su_C = MassConnector()
        self.su_H = MassConnector()
        self.ex_C = MassConnector()
        self.ex_H = MassConnector()

        self.Q_dot = HeatConnector()
        self.guesses = {}
        self.DT_pinch = -2
        
        self.DP_h = 0
        self.DP_c = 0
        
        self.h_hot = None
        self.h_cold = None
        self.T_hot = None
        self.T_cold = None

        self.eta_pinch = 1
        self.effectiveness = 1

    def get_required_inputs(self):
        self.sync_inputs()
        return ['P_su_H', 'T_su_H', 'm_dot_H', 'fluid_H', 'P_su_C', 'T_su_C', 'm_dot_C', 'fluid_C']
    
    def get_required_parameters(self):
        return ['eta_max', 'n_disc', 'Pinch_min']

    def pinch_residual(self, eta):
        self.counterflow_discretized(eta)
        return self.DT_pinch

    def find_Q_dot_max(self):
        # External enthalpy-based Q_dot_max
        h_h_id = _get_h_from_p_T(self.AS_H, self.su_H.fluid, self.su_H.p, self.su_C.T,
                                   su_connector=self.su_H)
        h_c_id = _get_h_from_p_T(self.AS_C, self.su_C.fluid, self.su_C.p, self.su_H.T,
                                   su_connector=self.su_C)

        Q_dot_maxh = self.su_H.m_dot * (self.su_H.h - h_h_id)
        Q_dot_maxc = self.su_C.m_dot * (h_c_id - self.su_C.h)

        self.Q_dot_max_ext = np.min([Q_dot_maxh, Q_dot_maxc])
        self.Q_dot_max = self.Q_dot_max_ext

        self.pinch_residual(1.0)
        if self.DT_pinch > 1e-3:
            self.eta_pinch = 1.0
            self.Q_dot_max_int = self.Q
            self.Q_dot_max = min(self.Q_dot_max_ext, self.Q_dot_max_int)
            return

        eta_low, eta_high = 0.0, 1.0
        tol_int = 1e-3

        for _ in range(30):
            eta_mid = 0.5 * (eta_low + eta_high)
            self.pinch_residual(eta_mid)

            if self.DT_pinch > tol_int:
                eta_low = eta_mid
            else:
                eta_high = eta_mid

            if abs(eta_high - eta_low) < 1e-3:
                break

        self.eta_pinch = eta_low
        self.pinch_residual(self.eta_pinch)
        self.Q_dot_max_int = self.Q
        self.Q_dot_max = min(self.Q_dot_max_ext, self.Q_dot_max_int)

    def counterflow_discretized(self, eta):
        n = self.params['n_disc']

        self.Q = eta * self.Q_dot_max
        Q_dot_seg = self.Q / n

        self.h_hot[0]  = self.su_H.h
        self.T_hot[0]  = self.su_H.T

        h_cold_out = self.su_C.h + self.Q / self.su_C.m_dot
        self.h_cold[0] = h_cold_out
        self.T_cold[0] = _get_T_from_h_p(self.AS_C, self.su_C.fluid,
                                           h_cold_out, self.su_C.p,
                                           su_connector=self.su_C)

        for i in range(n):
            self.h_hot[i+1] = self.h_hot[i] - Q_dot_seg / self.su_H.m_dot
            self.T_hot[i+1] = _get_T_from_h_p(self.AS_H, self.su_H.fluid,
                                                self.h_hot[i+1], self.p_hot[i+1],
                                                su_connector=self.su_H)

            self.h_cold[i+1] = self.h_cold[i] - Q_dot_seg / self.su_C.m_dot
            self.T_cold[i+1] = _get_T_from_h_p(self.AS_C, self.su_C.fluid,
                                                 self.h_cold[i+1], self.p_cold[i+1],
                                                 su_connector=self.su_C)

        self.DT_pinch = np.min(self.T_hot - self.T_cold)

    def solve(self):
        if self.su_H.m_dot is None:
            self.su_H.m_dot = self.su_C.m_dot
        if self.su_C.m_dot is None:
            self.su_C.m_dot = self.su_H.m_dot

        self.check_calculable()
        self.check_parametrized()

        if self.su_H.T < self.su_C.T:
            self.su_C, self.su_H = self.su_H, self.su_C

        self.DP_h = self.params.get('DP_h', 0)
        self.DP_c = self.params.get('DP_c', 0)

        if not self.calculable:
            print("HTX IS NOT CALCULABLE")
            return
        if not self.parametrized:
            print("HTX IS NOT PARAMETRIZED")
            return

        # Build AbstractState only for standard fluids; None for custom fluids
        self.AS_H = _get_abstract_state(self.su_H.fluid)
        self.AS_C = _get_abstract_state(self.su_C.fluid)

        if self.T_hot is None:
            n = self.params['n_disc']
            self.h_hot   = np.zeros(n + 1)
            self.h_cold  = np.zeros(n + 1)
            self.T_hot   = np.zeros(n + 1)
            self.T_cold  = np.zeros(n + 1)
            self.p_hot   = np.zeros(n + 1)
            self.p_cold  = np.zeros(n + 1)

            DP_h_disc = self.DP_h / n
            DP_c_disc = self.DP_c / n

            self.p_hot[0]  = self.su_H.p
            self.p_cold[0] = self.su_C.p

            for i in range(n):
                self.p_hot[i+1]  = self.p_hot[i]  - DP_h_disc
                self.p_cold[i+1] = self.p_cold[i] - DP_c_disc

        self.DT_pinch = -1
        self.find_Q_dot_max()
        self.DT_pinch = -1
        self.epsilon = self.params['eta_max']

        while self.DT_pinch <= self.params['Pinch_min']:
            self.counterflow_discretized(self.epsilon)

            if self.DT_pinch <= self.params['Pinch_min']:
                self.epsilon -= 0.01

                if self.epsilon <= 0:
                    self.solved = False
                    if self.print_flag:
                        print("No eta satisfies Pinch_min in HXEffCstDisc")
                    return

        self.update_connectors()
        self.solved = True

    def update_connectors(self):
        self.ex_C.set_fluid(self.su_C.fluid)
        self.ex_C.set_m_dot(self.su_C.m_dot)
        self.ex_C.set_h(self.su_C.h + self.Q / self.su_C.m_dot)
        self.ex_C.set_p(self.p_cold[-1])

        T_ex_C = _get_T_from_h_p(self.AS_C, self.su_C.fluid,
                                   self.ex_C.h, self.ex_C.p,
                                   su_connector=self.su_C)
        self.ex_C.T = T_ex_C

        # Hot side — bypass set_fluid for custom fluids
        h_ex_H = self.su_H.h - self.Q / self.su_H.m_dot
        if self.su_H.fluid in CUSTOM_FLUIDS:
            self.ex_H.fluid  = self.su_H.fluid
            self.ex_H.m_dot  = self.su_H.m_dot
            self.ex_H.p      = self.p_hot[-1]
            self.ex_H.h      = h_ex_H
            connector_class  = type(self.su_H)
            self.ex_H.T      = connector_class._T_from_h(h_ex_H)
            self.ex_H.cp     = connector_class._Cp(self.ex_H.T)
            self.ex_H.D      = connector_class._rho(self.ex_H.T)
            self.ex_H.state_known      = True
            self.ex_H.completely_known = True
        else:
            self.ex_H.set_fluid(self.su_H.fluid)
            self.ex_H.set_m_dot(self.su_H.m_dot)
            self.ex_H.set_h(h_ex_H)
            self.ex_H.set_p(self.p_hot[-1])
            T_ex_H = _get_T_from_h_p(self.AS_H, self.su_H.fluid,
                                       self.ex_H.h, self.ex_H.p,
                                       su_connector=self.su_H)
            self.ex_H.T = T_ex_H

        self.Q_dot.set_Q_dot(self.Q)

    def print_results(self):
        print("=== Heat Exchanger Results ===")
        print(f"Q: {self.Q_dot.Q_dot}")
        print("======================")

    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")
        print("======================")

    def plot_disc(self):
        import matplotlib.pyplot as plt
        x = np.arange(len(self.T_hot))
        plt.plot(x, self.T_hot,  label='Hot fluid [K]',  marker='o')
        plt.plot(x, self.T_cold, label='Cold fluid [K]', marker='s')
        plt.xlabel('Discretization segment')
        plt.ylabel('Temperature [K]')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()