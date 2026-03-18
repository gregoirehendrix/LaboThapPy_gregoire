# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:09:29 2026

@author: gregoire.hendrix
"""

from labothappy.machine.circuit_rec import RecursiveCircuit
from CoolProp.CoolProp import PropsSI

from labothappy.connector.mass_connector import MassConnector


from labothappy.component.compressor.compressor_csteff import CompressorCstEff 
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

def sCO2_basic(eta_cp, eta_ex, eta_heater, eta_cooler, HSource, P_low, P_high, T_c_su):
    sCO2 = RecursiveCircuit('CO2')
    
    #%% Create components
    compressor = CompressorCstEff()
    expander = ExpanderCstEff()
    gas_heater = HexCstEff()
    gas_cooler = HexCstEff()
    
    #%% Set components params
    compressor.set_parameters(eta_is = eta_cp)
    expander.set_parameters(eta_is = eta_ex)
    gas_heater.set_parameters(eta = eta_heater)
    gas_cooler.set_parameters(eta = eta_cooler)
    
    
    
    #%% Add components
    sCO2.add_component(compressor, "Compressor")
    sCO2.add_component(gas_heater, "Gas_Heater")
    sCO2.add_component(expander, "Expander")
    sCO2.add_component(gas_cooler, "Gas_Cooler")
    
    #%% Link components
    sCO2.link_components("Compressor", "m-ex", "Gas_Heater", "m-su_C")
    sCO2.link_components("Gas_Heater", "m-ex_C", "Expander", "m-su")
    sCO2.link_components("Expander", "m-ex", "Gas_Cooler", "m-su_H")
    sCO2.link_components("Gas_Cooler", "m-ex_H", "Compressor", "m-su")
    
    #%% Source & Sinks
    sCO2.add_source("Hot_air", HSource, sCO2.components["Gas_Heater"], "m-su_H")
    sCO2.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='Hot_air', P = HSource.p)
    
    sCO2.add_source("Cold_air", CSource, sCO2.components["Gas_Cooler"], "m-su_C")
    sCO2.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='Cold_air', P = CSource.p)
    
    #%%
    
    sCO2.set_cycle_guess(target='Compressor:su', m_dot = 15, T = T_c_su, p = P_low)
    sCO2.set_cycle_guess(target='Compressor:ex', p = P_high)
        
    sCO2.set_cycle_guess(target='Expander:su', p = P_high)
    sCO2.set_cycle_guess(target='Expander:ex', p = P_low)

    
    
    #%%
    sCO2.set_residual_variable(target='Expander:ex', variable='h', tolerance= 1e-3)
    sCO2.set_residual_variable(target='Expander:ex', variable='p', tolerance= 1e-3)

    sCO2.set_residual_variable(target='Compressor:ex', variable='h', tolerance= 1e-3)
    sCO2.set_residual_variable(target='Compressor:ex', variable='p', tolerance= 1e-3)

    sCO2.set_residual_variable(target='Gas_Heater:ex_C', variable='h', tolerance= 1e-3)
    sCO2.set_residual_variable(target='Gas_Heater:ex_C', variable='p', tolerance= 1e-3)

    sCO2.set_residual_variable(target='Gas_Cooler:ex_H', variable='h', tolerance= 1e-3)
    sCO2.set_residual_variable(target='Gas_Cooler:ex_H', variable='p', tolerance= 1e-3)
    

    return sCO2, compressor, expander, gas_heater, gas_cooler

#%%

if __name__ == "__main__":
    
    """--- Inputs ---"""
    
    # Pressure levels
    P_low_guess = 80*1e5            #Pa (=bar*1e5)
    P_high = 240*1e5                #Pa (=bar*1e5)
    
    # Temperature guess
    T_c_su_guess = 40 + 273.15
    # Isentropic efficiencies
    eta_heater = 0.968382                #effectiveness of hot heat exchanger
    eta_cooler = 0.9                #effectiveness of cold heat exchanger
    eta_ex = 0.9                    #isentropic efficiency of turbine/expander
    eta_cp = 0.8                    #isentropic efficiency of compressor
    
    # Hot source
    T_high = 565 + 273.15           #K
    p_high = 120*1e5                #Pa (=bar*1e5)
    m_dot_H = 20                     #kg/s
    
    # Cold source
    T_low = 25 + 273.15            #K
    p_low = 10*1e5                  #Pa (=bar*1e5)
    m_dot_C = 15                     #kg/s
    
    """--- Code (do not touch) ---"""
    
    HSource = MassConnector()
    HSource.set_properties(fluid = 'air', T = T_high, p = p_high, m_dot = m_dot_H)
    
    CSource = MassConnector()
    CSource.set_properties(fluid = 'air', T = T_low, p = p_low, m_dot = m_dot_C)
    
    sCO2, compressor, expander, gas_heater, gas_cooler = sCO2_basic(eta_cp, eta_ex, eta_heater, eta_cooler, HSource, P_low_guess, P_high, T_c_su_guess)
    sCO2.solve()
    
    compressor.print_work()
    expander.print_work()
    #gas_heater.print_results()
    #gas_cooler.print_results()
    
    
    
    sCO2.plot_cycle_Ts()
    
