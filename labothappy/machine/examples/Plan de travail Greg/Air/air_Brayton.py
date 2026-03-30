from labothappy.machine.circuit_it import IterativeCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff
from labothappy.connector.solar_salt_connector import SolarSaltConnector

# --- Parameters ---
fluid       = 'Air'
T_amb       = 13.2 + 273.15
P_atm       = 1.01325e5
PR          = 11
eta_comp    = 0.8
eta_turb    = 0.9
eta_HX      = 1
P_salt      = 1e5
T_salt_cold = 290 + 273.15

T_salt_limit = 565 + 273.15   # K — switch threshold
T_hot_su     = 850 + 273.15   # K — change this to go above salt limit

# --- Hot source selection ---
if T_hot_su <= T_salt_limit:
    hot_source = SolarSaltConnector()
    hot_source.set_properties(T=T_hot_su, p=P_salt, m_dot=10.0)
    hot_fluid  = 'SolarSalt'
    #print(f"Hot source : Solar Salt at {T_hot_su - 273.15:.1f} °C")
else:
    hot_source = MassConnector()
    hot_source.set_properties(fluid='Air', T=T_hot_su, p=P_salt, m_dot=10.0)
    hot_fluid  = 'Air'
    #print(f"Hot source : Hot Air at {T_hot_su - 273.15:.1f} °C")

def brayton_simple():
    # --- Circuit ---
    cycle = IterativeCircuit(fluid=fluid)
    
    Compressor = CompressorCstEff()
    Heater     = HexCstEff()
    Turbine    = ExpanderCstEff()
    
    Compressor.set_parameters(eta_is=eta_comp)
    Turbine.set_parameters(eta_is=eta_turb)
    Heater.set_parameters(eta=eta_HX)
    
    cycle.add_component(Compressor, "Compressor")
    cycle.add_component(Heater,     "Heater")
    cycle.add_component(Turbine,    "Turbine")
    
    # --- Links (air side) ---
    cycle.link_components("Compressor", "m-ex",   "Heater",  "m-su_C")
    cycle.link_components("Heater",     "m-ex_C", "Turbine", "m-su")
    
    # --- Air inlet source ---
    air_in = MassConnector()
    cycle.add_source("AirInlet", air_in, cycle.components["Compressor"], "m-su")
    cycle.set_source_properties(
        target = "AirInlet",
        fluid  = fluid,
        T      = T_amb,
        P      = P_atm,
        m_dot  = 1.0
    )
    
    # --- Hot source ---
    cycle.add_source("HotSource", hot_source, cycle.components["Heater"], "m-su_H")
    cycle.set_source_properties(
        target = "HotSource",
        fluid  = hot_fluid,
        T      = T_hot_su,
        P      = P_salt,
        m_dot  = 10.0
    )
    
    # --- Imposed pressures ---
    cycle.set_cycle_input(target="Compressor:ex", p=P_atm * PR)
    cycle.set_cycle_input(target="Turbine:ex",    p=P_atm)
    
    # --- Sequential solve ---
    cycle._build_solve_order()
    for name in cycle.solve_start_components:
        cycle.components[name].solve()



# --- Results ---
try:
    W_comp = Compressor.W.W_dot
except NameError:
    print("J'étais en train de réarranger par fonctions")
    sys.exit()
W_turb = Turbine.W.W_dot
W_net  = W_turb - W_comp
Q_in   = Heater.Q.Q_dot
eta_th = W_net / Q_in

# --- m_dot_hot source from energy balance ---
if T_hot_su <= T_salt_limit:
    Cp_salt_avg = SolarSaltConnector._Cp((T_hot_su + T_salt_cold) / 2)
    m_dot_hot   = Q_in / (Cp_salt_avg * (T_hot_su - T_salt_cold))
    hot_label   = "m_dot_salt"
else:
    from CoolProp.CoolProp import PropsSI
    T_hot_ex    = Heater.ex_H.T
    h_hot_in    = PropsSI('H', 'T', T_hot_su,  'P', P_salt, 'Air')
    h_hot_ex    = PropsSI('H', 'T', T_hot_ex,  'P', P_salt, 'Air')
    m_dot_hot   = Q_in / (h_hot_in - h_hot_ex)
    hot_label   = "m_dot_air_hot"

print(f"{'─'*45}")
print(f"  Simple Open Air Brayton Cycle")
print(f"  Hot source : {hot_fluid} at {T_hot_su - 273.15:.1f} °C")
print(f"{'─'*45}")
print(f"  PR             = {PR:.1f}  [-]")
print(f"  TIT            = {Turbine.su.T - 273.15:.1f}  [°C]")
print(f"  T_ex_comp      = {Compressor.ex.T - 273.15:.1f}  [°C]")
print(f"  T_ex_turb      = {Turbine.ex.T - 273.15:.1f}  [°C]")
print(f"  dT_turb        = {Turbine.su.T - Turbine.ex.T:.1f}  [°C]")
print(f"  T_ex_heater_H  = {Heater.ex_H.T - 273.15:.1f}  [°C]")
print(f"  W_comp         = {W_comp/1e3:.2f}  [kW per kg/s air]")
print(f"  W_turb         = {W_turb/1e3:.2f}  [kW per kg/s air]")
print(f"  W_net          = {W_net/1e3:.2f}  [kW per kg/s air]")
print(f"  Q_in           = {Q_in/1e3:.2f}  [kW per kg/s air]")
print(f"  eta_th         = {eta_th*100:.1f}  [%]")
print(f"  {hot_label:<15}= {m_dot_hot:.4f}  [kg/s per kg/s air]")
print(f"{'─'*45}")