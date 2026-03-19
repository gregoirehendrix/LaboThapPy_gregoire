# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 08:42:09 2026

@author: gregoire.hendrix
"""

from labothappy.machine.circuit_rec import RecursiveCircuit

from labothappy.connector.mass_connector import MassConnector

from labothappy.component.compressor.compressor_csteff import CompressorCstEff 
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

def air_basic(eta_cp, eta_tb, eta_hx, HSource, AirAmb, m_dot_air, P_high):
    try:
        
        air = RecursiveCircuit('air')
        
        #%% Ceate components
        compressor = CompressorCstEff()
        turbine = ExpanderCstEff()
        heater = HexCstEff()
        
        #%% Set components params
        compressor.set_parameters(eta_is = eta_cp)
        turbine.set_parameters(eta_is = eta_tb)
        heater.set_parameters(eta = eta_hx)
        
        #%% Add components
        air.add_component(compressor, 'Compressor')
        air.add_component(turbine, 'Turbine')
        air.add_component(heater, 'Heater')
        
        #%% Link components
        air.link_components('Compressor', 'm-ex', 'Heater', 'm-su_C')
        air.link_components('Heater', 'm-ex_C', 'Turbine', 'm-su')
        
        #%% Add sources & sinks
        air.add_source('HotAir', HSource, air.components["Heater"], 'm-su_H')
        air.add_source('AmbientAir', AirAmb, air.components["Compressor"], 'm-su')
        air.add_source('AmbientAir', AirAmb, air.components["Turbine"], 'm-ex')
        
        #%% Guesses
        #Compressor Guesses        
        air.set_cycle_guess(target='Compressor:su', m_dot = m_dot_air)
        air.set_cycle_guess(target='Compressor:ex', p = P_high)
        # Turbine guesses    
        air.set_cycle_guess(target='Turbine:su', p = P_high)
        
        #%% Tolerance limit
        air.set_residual_variable(target='Compressor:ex', variable='h', tolerance= 1e-3)
        air.set_residual_variable(target='Compressor:ex', variable='p', tolerance= 1e-3)
        air.set_residual_variable(target='Heater:ex_C', variable='h', tolerance= 1e-3)
        air.set_residual_variable(target='Heater:ex_C', variable='p', tolerance= 1e-3)
        air.set_residual_variable(target='Turbine:ex', variable='h', tolerance= 1e-3)
        air.set_residual_variable(target='Turbine:ex', variable='p', tolerance= 1e-3)     
    
    
    except Exception as e:
            print("Une erreur est survenue :", e)
            raise

    else:
        print("Le code à tourné sans erreur !")

    
    
    return air

#%%

if __name__ == "__main__":
    
    """ --- Inputs --- """
    eta_cp = 0.8
    eta_tb = 0.9
    eta_hx = 0.9
    
    # Guesses
    m_dot_air = 10
    P_high = 10*1e5
    
    # Hot source
    T_high = 700 + 273.15
    p_high = 10*1e5
    m_dot_H = 10
    
    # Ambient air
    T_amb = 25 + 273.15
    p_amb = 101325
    
    HSource = MassConnector()
    HSource.set_properties(fluid = 'air', T = T_high, P = p_high, m_dot = m_dot_H)
    
    AirAmb = MassConnector()
    AirAmb.set_properties(fluid = 'air', T = T_amb, P = p_amb)
    
    air_Brayton = air_basic(eta_cp, eta_tb, eta_hx, HSource, AirAmb, m_dot_air, P_high)
    air_Brayton.solve()
    