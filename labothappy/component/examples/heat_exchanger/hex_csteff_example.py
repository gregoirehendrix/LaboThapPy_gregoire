from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

"Simple test - Model for recuperators in simple thermodynamic studies "

HTX = HexCstEff()

# Set input conditions
HTX.set_inputs(
    fluid_C='air',
    T_su_C=273.15 + 100,
    m_dot_C=1,
    P_su_C=10e5,

    fluid_H='CO2',
    T_su_H=273.15 + 400,
    m_dot_H=1,
    P_su_H=10e5,
)

# Set heat exchanger parameters
HTX.set_parameters(**{
    'eta': 0.95,
})

# Solve the model
HTX.solve()

# Output results
HTX.print_results()
HTX.print_states_connectors()
