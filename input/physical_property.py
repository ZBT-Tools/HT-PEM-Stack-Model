""" This file contains the physical settings"""

"""Electrochemistry Settings"""
# cathode volumetric exchange current density [A/m^3]
cathode_exchange_current_density = 2.3e3
# anode volumetric exchange current density [A/m^3]
anode_exchange_current_density = 0.817e9
# cathode layer diffusion coefficient [m^2/s]
cathode_layer_diffusion_coefficient = 1.36e-8
# anode layer diffusion coefficient [m^2/s]
anode_layer_diffusion_coefficient = 1.36e-8
# cathode gdl diffusion coefficient [m^2/s]
cathode_gdl_diffusion_coefficient = 2.59e-6
# anode gdl diffusion coefficient [m^2/s]
anode_gdl_diffusion_coefficient = 2.59e-6
# cathode layer proton conductivity [Ohm^-1m^-1]
cathode_proton_conductivity = 3.e0
# anode layer proton conductivity [Ohm^-1m^-1]
anode_proton_conductivity = 3.e0
# cathode tafel slope [V]
cathode_tafel_slope = 0.03
# anode tafel slope [V]
anode_tafel_slope = 0.03
# thermo neutral open circuit voltage [V]
v_thermo_neutral = 1.28
# fitted vapour mass transport coefficient [m/s]
vapour_mass_transport_coefficient = 0.62e-5
# molar concentration of membrane acid groups [mol/m³]
molar_membrane_acid_group_concentration = 1.2e3


"""Electric Settings"""
# membrane basic resistance [Ohm/m^2]  name ???
membrane_basic_resistance = 0.33
# membrane temperature slope
membrane_temperature_resistance = 0.0007
# bipolar plate resistivity [Ohm/m]
bipolar_plate_resistivity = 2.e-6


"""Thermal Settings"""
# thermal conductivity of the bipolar plate through plane  [W/(mK)]
k_p = 1.e2
# thermal conductivity of the bipolar plate in plane  [W/(mK)]
ki_p = 1.e2
# thermal conductivity of the gde through plane  [W/(mK)]
k_g = 1.e0
# thermal conductivity of the gde in plane  [W/(mK)]
ki_g = 1.e0
# thermal conductivity of the membrane through plane  [W/(mK)]
k_m = .26e0
# thermal conductivity of the membrane in plane  [W/(mK)]
ki_m = .26e0
# coolant heat capacity [J/(kgK)]
cp_coolant = 4.e3
# convection coefficient of the coolant/channel [W/(Km^2)]
convection_coefficient_coolant = 4.e3
# convection coefficient of the stack walls/environment [W/(Km^2)]
convection_coefficient_stack = 5.e0
# latent heat of vaporization [J/mol]
enthalpy_vaporization = 45.4e3

"""Fluid Mechanic Settings"""
# geometrical pressure loss coefficient of the manifold header
manifold_pressure_loss_coefficient = 5.e5
# bend pressure loss coefficient of the channel bends
bend_pressure_loss_coefficient = 0.1
# number of channel bends
channel_bends = 48
# cathode channel gas flow direction
cathode_channel_flow_direction = True
# anode channel gas flow direction
anode_channel_flow_direction = False

"""Humidification"""
# cathode inlet gas relative humidity
cathode_inlet_humidity = 0.
# anode inlet gas relative humidity
anode_inlet_humidity = 0.

"""Species settings"""
# oxygen concentration at the cathode inlet
oxygen_inlet_concentration = 0.21
# hydrogen concentration at the anode inlet
hydrogen_inlet_concentration = 0.5

