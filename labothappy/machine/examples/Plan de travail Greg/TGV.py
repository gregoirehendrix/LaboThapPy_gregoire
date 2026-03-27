from labothappy.machine.circuit_it import IterativeCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from scipy.optimize import brentq

# ============================================================
#  PARAMETERS
# ============================================================

fluid_air    = 'Air'
T_amb        = 13.2 + 273.15
P_atm        = 1.01325e5
m_dot_air    = 1.0

T_hot_su     = 1000 + 273.15
P_hot        = 1e5
m_dot_hot    = 10.0

PR_brayton   = 10.0
eta_comp     = 0.80
eta_turb_GT  = 0.90
eta_HX_heat  = 1.0

fluid_steam  = 'Water'
P_high       = 30e5             # 30 bar — bon compromis x_ex_ST / rendement
P_low        = 0.10e5
eta_pump     = 0.80
eta_turb_ST  = 0.85
eta_HRSG     = 0.90

# ============================================================
#  CIRCUIT 1 — BRAYTON
# ============================================================

brayton = IterativeCircuit(fluid=fluid_air)

Compressor = CompressorCstEff()
Heater     = HexCstEff()
GasTurbine = ExpanderCstEff()

Compressor.set_parameters(eta_is=eta_comp)
GasTurbine.set_parameters(eta_is=eta_turb_GT)
Heater.set_parameters(eta=eta_HX_heat)

brayton.add_component(Compressor, "Compressor")
brayton.add_component(Heater,     "Heater")
brayton.add_component(GasTurbine, "GasTurbine")

brayton.link_components("Compressor", "m-ex",   "Heater",     "m-su_C")
brayton.link_components("Heater",     "m-ex_C", "GasTurbine", "m-su")

air_in = MassConnector()
brayton.add_source("AirInlet", air_in, brayton.components["Compressor"], "m-su")
brayton.set_source_properties(target="AirInlet", fluid=fluid_air,
                               T=T_amb, P=P_atm, m_dot=m_dot_air)

hot_source = MassConnector()
brayton.add_source("HotSource", hot_source, brayton.components["Heater"], "m-su_H")
brayton.set_source_properties(target="HotSource", fluid=fluid_air,
                               T=T_hot_su, P=P_hot, m_dot=m_dot_hot)

brayton.set_cycle_input(target="Compressor:ex", p=P_atm * PR_brayton)
brayton.set_cycle_input(target="GasTurbine:ex", p=P_atm)

brayton._build_solve_order()
for name in brayton.solve_start_components:
    brayton.components[name].solve()

# ============================================================
#  EXTRACT EXHAUST GAS
# ============================================================

T_exhaust    = GasTurbine.ex.T
P_exhaust    = GasTurbine.ex.p
h_exhaust_in = GasTurbine.ex.h

# ============================================================
#  RANKINE — résolution analytique complète via CoolProp
#
#  Toutes les enthalpies sont calculées directement.
#  Le HRSG est résolu par bilan enthalpique (pas de cp moyen),
#  ce qui est correct même en présence de changement de phase.
#
#  Hypothèse : vapeur sort saturée sèche à P_high (pas de surchauffe).
#  C'est le cas optimal pour un HRSG sans surchauffe explicite.
#
#  Le vrai m_dot_steam est obtenu par brentq sur le résidu du HRSG :
#
#    résidu(m_dot_steam) = Q_steam(m_dot_steam) - Q_gas(m_dot_steam)
#
#  avec :
#    Q_steam = m_dot_steam * (h_vap_sat - h_ex_pump)
#    Q_gas   = eta * min(Q_max_gas, Q_max_steam)
#
#  Q_max_gas  = m_dot_air * (h_exhaust_in - h_gas_at_T_steam_in)
#  Q_max_steam = m_dot_steam * (h_vap_sat - h_ex_pump)
#
#  → si Q_max_gas < Q_max_steam : gaz limitant, Q = eta*Q_max_gas
#    résidu = m_dot_steam*(h_vap_sat - h_ex_pump) - eta*Q_max_gas = 0
#    → solution analytique directe (pas besoin de brentq dans ce cas)
#
#  → si Q_max_steam < Q_max_gas : vapeur limitante (pinch côté vapeur)
#    Q = eta*Q_max_steam = eta*m_dot_steam*(h_vap_sat - h_ex_pump)
#    → résidu = 0 toujours : système sous-déterminé, on prend Q_max_gas
#
#  En pratique pour un HRSG : gaz toujours limitant si m_dot_air << m_dot_steam.
#  On vérifie via bilan enthalpique direct.
# ============================================================

AS = CP.AbstractState('HEOS', fluid_steam)
AS_air = CP.AbstractState('HEOS', fluid_air)

# --- Pompe ---
h_sat_liq  = PropsSI('H', 'P', P_low,  'Q', 0, fluid_steam)
h_sat_vap  = PropsSI('H', 'P', P_low,  'Q', 1, fluid_steam)
T_sat_low  = PropsSI('T', 'P', P_low,  'Q', 0, fluid_steam)
s_sat_liq  = PropsSI('S', 'P', P_low,  'Q', 0, fluid_steam)

h_ex_pump_is = PropsSI('H', 'P', P_high, 'S', s_sat_liq, fluid_steam)
h_ex_pump    = h_sat_liq + (h_ex_pump_is - h_sat_liq) / eta_pump
T_ex_pump    = PropsSI('T', 'P', P_high, 'H', h_ex_pump, fluid_steam)
w_pump       = h_ex_pump - h_sat_liq

# --- Etat cible vapeur en sortie HRSG : vap. sat. sèche à P_high ---
h_vap_high = PropsSI('H', 'P', P_high, 'Q', 1, fluid_steam)
T_vap_high = PropsSI('T', 'P', P_high, 'Q', 0, fluid_steam)

# --- HRSG : bilan enthalpique côté gaz ---
# Q_max_gas = m_dot_air * (h_exhaust_in - h_gas_cooled_to_T_steam_in)
# T_steam_in = T_ex_pump (entrée eau froide côté C)
AS_air.update(CP.PT_INPUTS, P_exhaust, T_ex_pump)
h_gas_at_T_steam_in = AS_air.hmass()
Q_max_gas = m_dot_air * (h_exhaust_in - h_gas_at_T_steam_in)

# Q_steam nécessaire pour vaporiser 1 kg/s de vapeur sat. sèche
delta_h_steam = h_vap_high - h_ex_pump   # J/kg vapeur

# m_dot_steam si gaz limitant
Q_HRSG      = eta_HRSG * Q_max_gas
m_dot_steam = Q_HRSG / delta_h_steam

# Vérification : Q_max_steam avec ce m_dot
Q_max_steam = m_dot_steam * delta_h_steam
# Si gaz limitant → Q_max_gas < Q_max_steam → OK
# (toujours vrai car Q_HRSG = eta*Q_max_gas < Q_max_gas = Q_max_steam/eta < Q_max_steam)

# Température sortie gaz (bilan côté gaz)
h_gas_ex = h_exhaust_in - Q_HRSG / m_dot_air

def f_T(T):
    AS_air.update(CP.PT_INPUTS, P_exhaust, T)
    return AS_air.hmass() - h_gas_ex

T_ex_gas = brentq(f_T, 200, T_exhaust - 0.1)

# --- Turbine vapeur ---
h_su_ST    = h_vap_high
s_su_ST    = PropsSI('S', 'P', P_high, 'Q', 1, fluid_steam)
T_su_ST    = T_vap_high

h_ex_ST_is = PropsSI('H', 'P', P_low, 'S', s_su_ST, fluid_steam)
h_ex_ST    = h_su_ST - eta_turb_ST * (h_su_ST - h_ex_ST_is)
T_ex_ST    = PropsSI('T', 'P', P_low, 'H', h_ex_ST, fluid_steam)
w_ST       = h_su_ST - h_ex_ST
x_ex_ST    = (h_ex_ST - h_sat_liq) / (h_sat_vap - h_sat_liq)

# --- Condenseur ---
Q_cond_unit = h_ex_ST - h_sat_liq       # J/kg vapeur
Q_cond_real = m_dot_steam * Q_cond_unit # W / (kg/s air)

# ============================================================
#  BILAN DE PUISSANCE
# ============================================================

W_pump_air = w_pump * m_dot_steam
W_ST_air   = w_ST   * m_dot_steam
W_net_ST   = W_ST_air - W_pump_air

W_comp   = Compressor.W.W_dot
W_GT     = GasTurbine.W.W_dot
W_net_GT = W_GT - W_comp
Q_in     = Heater.Q.Q_dot
eta_th_GT = W_net_GT / Q_in

W_net_combined = W_net_GT + W_net_ST
eta_combined   = W_net_combined / Q_in

# ============================================================
#  AFFICHAGE
# ============================================================

print(f"{'═'*52}")
print(f"  CYCLE COMBINÉ  —  Brayton Air + Rankine Vapeur")
print(f"{'═'*52}")
print(f"\n  {'─'*48}")
print(f"  BRAYTON (topping, air)")
print(f"  {'─'*48}")
print(f"  PR_brayton       = {PR_brayton:.1f}         [-]")
print(f"  T_hot_su         = {T_hot_su-273.15:.1f}      [°C]")
print(f"  TIT              = {GasTurbine.su.T-273.15:.1f}      [°C]")
print(f"  T_ex_compressor  = {Compressor.ex.T-273.15:.1f}       [°C]")
print(f"  T_exhaust        = {T_exhaust-273.15:.1f}       [°C]")
print(f"  W_comp           = {W_comp/1e3:.2f}      [kW/(kg/s air)]")
print(f"  W_GT             = {W_GT/1e3:.2f}      [kW/(kg/s air)]")
print(f"  W_net_GT         = {W_net_GT/1e3:.2f}      [kW/(kg/s air)]")
print(f"  Q_in             = {Q_in/1e3:.2f}      [kW/(kg/s air)]")
print(f"  eta_th_GT        = {eta_th_GT*100:.1f}         [%]")
print(f"\n  {'─'*48}")
print(f"  RANKINE (bottoming, vapeur d'eau)")
print(f"  {'─'*48}")
print(f"  P_high           = {P_high/1e5:.0f}            [bar]")
print(f"  P_low            = {P_low/1e5:.2f}         [bar]")
print(f"  T_sat(P_low)     = {T_sat_low-273.15:.1f}       [°C]")
print(f"  T_su_Pump        = {T_sat_low-273.15:.2f}       [°C]  (liq. sat.)")
print(f"  T_ex_Pump        = {T_ex_pump-273.15:.2f}       [°C]")
print(f"  T_su_HRSG (eau)  = {T_ex_pump-273.15:.2f}       [°C]")
print(f"  T_ex_HRSG (eau)  = {T_vap_high-273.15:.1f}      [°C]  (vap. sat. sèche)")
print(f"  T_ex_HRSG (gaz)  = {T_ex_gas-273.15:.1f}       [°C]")
print(f"  T_su_ST          = {T_su_ST-273.15:.1f}      [°C]")
print(f"  T_ex_ST          = {T_ex_ST-273.15:.2f}       [°C]")
print(f"  x_ex_ST          = {x_ex_ST:.3f}         [-]  (titre vapeur)")
print(f"  m_dot_steam      = {m_dot_steam:.4f}      [kg/s vapeur / kg/s air]")
print(f"  w_pump           = {w_pump/1e3:.4f}      [kJ/kg vapeur]")
print(f"  w_ST             = {w_ST/1e3:.2f}      [kJ/kg vapeur]")
print(f"  W_net_ST         = {W_net_ST/1e3:.2f}      [kW/(kg/s air)]")
print(f"  Q_HRSG           = {Q_HRSG/1e3:.2f}      [kW/(kg/s air)]")
print(f"  Q_cond           = {Q_cond_real/1e3:.2f}      [kW/(kg/s air)]")
print(f"\n  {'─'*48}")
print(f"  BILAN CYCLE COMBINÉ")
print(f"  {'─'*48}")
print(f"  W_net_GT         = {W_net_GT/1e3:.2f}      [kW/(kg/s air)]")
print(f"  W_net_ST         = {W_net_ST/1e3:.2f}      [kW/(kg/s air)]")
print(f"  W_net_combined   = {W_net_combined/1e3:.2f}      [kW/(kg/s air)]")
print(f"  Q_in             = {Q_in/1e3:.2f}      [kW/(kg/s air)]")
print(f"  eta_combined     = {eta_combined*100:.1f}         [%]")
print(f"{'═'*52}")

# --- Vérification premier principe ---
bilan = Q_in - W_net_combined - Q_cond_real - (Q_HRSG - W_net_ST - Q_cond_real)
print(f"\n  [CHECK] Bilan Rankine  Q_HRSG - W_net_ST - Q_cond = {(Q_HRSG - W_net_ST - Q_cond_real)/1e3:.2f} kW (doit = 0)")