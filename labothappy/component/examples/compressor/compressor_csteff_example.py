from labothappy.component.compressor.compressor_csteff import CompressorCstEff

# Example usage
CP = CompressorCstEff()
# CP.print_setup()

"If the inputs are not set directly BUT through the connectors"
# CP.su.set_properties(P=319296.5575177148, T=331.033964665788, fluid='R1233ZDE', m_dot = 0.1)
# CP.ex.set_properties(P=606240.1433176235)

"If the inputs are set directly"
CP.set_inputs(
    P_su=80*1e5,
    T_su=273.15+50,
    P_ex=240*1e5,
    fluid='CO2',  # Make sure to include fluid information
    m_dot=100  # Mass flow rate
)
CP.set_parameters(eta_is=0.8)
# CP.print_setup()



CP.solve()
#CP.print_results()
CP.print_work()

fig = CP.plot_Ts()
fig.show()
