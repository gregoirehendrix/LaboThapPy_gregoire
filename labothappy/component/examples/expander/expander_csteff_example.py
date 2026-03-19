from labothappy.component.expander.expander_csteff import ExpanderCstEff

# Example usage
EXP = ExpanderCstEff()
EXP.print_setup()

set_prop = False

# "If the inputs are not set directly BUT throught the connectors"
if set_prop :
    EXP.su.set_properties(
        P=955214.9, 
        T=374.18, 
        fluid='R134a', 
        m_dot = 0.1)
    EXP.ex.set_properties(P=293940.1)
else:
    
    EXP.set_inputs(
        P_su=240*1e5,
        T_su=273.15+550,
        P_ex=80*1e5,
        fluid='CO2',  # Make sure to include fluid information
        m_dot=20  # Mass flow rate
    )

EXP.set_parameters(eta_is=0.9)
EXP.print_setup()

EXP.solve()
EXP.print_results()
EXP.print_work()

fig = EXP.plot_Ts()
fig.show()

