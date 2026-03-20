from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

import CoolProp.CoolProp as CP

# Fluids not available in CoolProp — properties must come from the connector directly
CUSTOM_FLUIDS = {'SolarSalt'}

class HexCstEff(BaseComponent):
    """
    **Component**: Counterflow Heat Exchanger with Constant Effectiveness (HXEffCst)
    
    **Model**: Simplified Heat Exchanger Model with Constant Effectiveness
    
    **Description**:
    
        This model simulates a counterflow heat exchanger with a fixed effectiveness (η). It estimates the maximum possible heat transfer based on inlet conditions and applies the given effectiveness to compute the actual heat transfer rate. The model assumes a simplified configuration without discretization or detailed internal states. It is suitable for quick, steady-state evaluations in system-level models.
    
    **Assumptions**:
    
        - Steady-state operation.
        - Constant effectiveness (η) is applied across the entire heat exchanger.
        - No phase change or pressure drop effects are modeled.
        - Uniform flow distribution on both hot and cold sides.
        - Fluid properties are retrieved from CoolProp.
        - For custom fluids (e.g. SolarSalt), properties are taken from the connector directly.
    
    **Connectors**:
    
        su_C (MassConnector): Mass connector for the cold-side supply.
        
        su_H (MassConnector): Mass connector for the hot-side supply.
        
        ex_C (MassConnector): Mass connector for the cold-side exhaust.
        
        ex_H (MassConnector): Mass connector for the hot-side exhaust.
        
        Q (HeatConnector): Heat transfer connector representing the total exchanged heat.
    
    **Parameters**:
    
        eta (float): Effectiveness of the heat exchanger [-].
    
    **Inputs**:
    
        fluid_H (str): Hot-side fluid.
        
        h_su_H (float): Hot-side inlet specific enthalpy [J/kg].
        
        P_su_H (float): Hot-side inlet pressure [Pa].
        
        m_dot_H (float): Hot-side mass flow rate [kg/s].
    
        fluid_C (str): Cold-side fluid.
        
        h_su_C (float): Cold-side inlet specific enthalpy [J/kg].
        
        P_su_C (float): Cold-side inlet pressure [Pa].
        
        m_dot_C (float): Cold-side mass flow rate [kg/s].
    
    **Outputs**:
    
        h_ex_C: Cold-Side Exhaust specific enthalpy at outlet [J/kg].
        
        P_ex_C: Cold-Side Exhaust pressure at outlet [Pa].
                    
        h_ex_H: Hot-Side specific enthalpy at outlet [J/kg].
        
        P_ex_H: Hot-Side pressure at outlet [Pa].
        
        Q_dot: Total heat transfer rate across the exchanger [W].

    """

    def __init__(self):
        super().__init__()
        self.su_C = MassConnector()
        self.su_H = MassConnector()
        self.ex_C = MassConnector()
        self.ex_H = MassConnector()

        self.Q = HeatConnector()
        self.guesses = {}

    def get_required_inputs(self):
        return ['fluid_C', 'h_su_C', 'P_su_C', 'm_dot_C', 'fluid_H', 'h_su_H', 'P_su_H', 'm_dot_H']
    
    def get_required_parameters(self):
        return ['eta']

    def solve(self):

        if self.su_H.m_dot is None: 
            self.su_H.m_dot = self.su_C.m_dot
        if self.su_C.m_dot is None:
            self.su_C.m_dot = self.su_H.m_dot    

        self.check_calculable()
        self.check_parametrized()
        
        self.DP_c = self.params.get('DP_c', 0)
        self.DP_h = self.params.get('DP_h', 0)

        # --- Path 1: enthalpy-based Q_dot_max via CoolProp (standard fluids, most accurate) ---
        try:
            self.AS_C = CP.AbstractState('HEOS', self.su_C.fluid)
            self.AS_H = CP.AbstractState('HEOS', self.su_H.fluid)

            self.AS_H.update(CP.PT_INPUTS, self.su_H.p, self.su_C.T)
            H_h_id = self.AS_H.hmass()
            self.AS_C.update(CP.PT_INPUTS, self.su_C.p, self.su_H.T)
            H_c_id = self.AS_C.hmass()

            Q_dot_maxh = self.su_H.m_dot * abs(H_h_id - self.su_H.h)
            Q_dot_maxc = self.su_C.m_dot * abs(H_c_id - self.su_C.h)
            Q_dot_max  = min(Q_dot_maxh, Q_dot_maxc)

            Q_dot = self.params['eta'] * Q_dot_max
            self.update_connectors(Q_dot)
            self.solved = True
            return

        except Exception:
            pass

        # --- Path 2: Cp-based Q_dot_max via CoolProp (fallback for standard fluids) ---
        try:
            self.AS_C = CP.AbstractState('HEOS', self.su_C.fluid)
            self.AS_H = CP.AbstractState('HEOS', self.su_H.fluid)

            self.AS_C.update(CP.PT_INPUTS, self.su_C.p, self.su_C.T)
            cp_c = self.AS_C.cpmass()
            self.AS_H.update(CP.PT_INPUTS, self.su_H.p, self.su_H.T)
            cp_h = self.AS_H.cpmass()

            C_dot_min = min(self.su_C.m_dot * cp_c, self.su_H.m_dot * cp_h)
            Q_dot_max = C_dot_min * (self.su_H.T - self.su_C.T)

            Q_dot = self.params['eta'] * Q_dot_max
            self.update_connectors(Q_dot)
            self.solved = True
            return

        except Exception:
            pass

        # --- Path 3: Cp from connector (custom fluids — e.g. SolarSalt) ---
        # Requires su_H.cp and su_C.cp to be pre-computed by the custom connector.
        cp_c = self.su_C.cp
        cp_h = self.su_H.cp

        C_dot_min = min(self.su_C.m_dot * cp_c, self.su_H.m_dot * cp_h)
        Q_dot_max = C_dot_min * (self.su_H.T - self.su_C.T)

        Q_dot = self.params['eta'] * Q_dot_max
        self.update_connectors(Q_dot)
        self.solved = True

    def update_connectors(self, Q_dot):

        # --- Cold side ---
        self.ex_C.set_fluid(self.su_C.fluid)
        self.ex_C.set_m_dot(self.su_C.m_dot)
        self.ex_C.set_h(self.su_C.h + Q_dot / self.su_C.m_dot)
        self.ex_C.set_p(self.su_C.p - self.DP_c)

        # --- Hot side ---
        # For custom fluids (not in CoolProp), bypass set_fluid to avoid CoolProp lookup.
        # Properties are recovered from the source connector's own correlations.
        if self.su_H.fluid in CUSTOM_FLUIDS:
            self.ex_H.fluid  = self.su_H.fluid
            self.ex_H.m_dot  = self.su_H.m_dot
            self.ex_H.p      = self.su_H.p - self.DP_h
            h_ex_H           = self.su_H.h - Q_dot / self.su_H.m_dot
            self.ex_H.h      = h_ex_H
            # Recover T, cp, rho from outlet enthalpy using the source connector's methods
            self.ex_H.T      = type(self.su_H)._T_from_h(h_ex_H)
            self.ex_H.cp     = type(self.su_H)._Cp(self.ex_H.T)
            self.ex_H.D      = type(self.su_H)._rho(self.ex_H.T)
            self.ex_H.state_known     = True
            self.ex_H.completely_known = True
        else:
            self.ex_H.set_fluid(self.su_H.fluid)
            self.ex_H.set_m_dot(self.su_H.m_dot)
            self.ex_H.set_h(self.su_H.h - Q_dot / self.su_H.m_dot)
            self.ex_H.set_p(self.su_H.p - self.DP_h)

        # --- Heat connector ---
        self.Q.set_Q_dot(Q_dot)

    def print_results(self):
        print("=== Heat Exchanger Results ===")
        print(f"Q: {self.Q.Q_dot/1000:.3f} kW")
        print("======================")

    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        print("Connectors:")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q.Q_dot}")
        print("======================")