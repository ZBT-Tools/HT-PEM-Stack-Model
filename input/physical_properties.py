""" This file contains the physical settings"""

"""Electrochemistry Settings"""
# volumetric exchange current density of the cathode[A/m^3]
exchange_current_density_cathode = 1.4e5
# volumetric exchange current density of the anode [A/m^3]
exchange_current_density_anode = 0.817e9
# oxygen catalyst layer diffusion coefficient [m^2/s]
oxygen_catalyst_layer_diffusion_coefficient = 1.36e-8
# hydrogen catalyst layer diffusion coefficient [m^2/s]
hydrogen_catalyst_layer_diffusion_coefficient = 5.e-8
# oxygen gas diffusion layer diffusion coefficient [[m^2/s]
oxygen_gas_diffusion_layer_diffusion_coefficient = 2.59e-6
# hydrogen gas diffusion layer diffusion coefficient [m^2/s]
hydrogen_diffusion_layer_diffusion_coefficient = 9.52e-6
# catalyst layer proton conductivity of the cathode [Ohm^-1/m]
catalyst_layer_proton_conductivity_cathode = 3.e0
# catalyst layer proton conductivity of the anode [Ohm^-1/m]
catalyst_layer_proton_conductivity_anode = 3.e0
# tafel slope of the cathode [V]
tafel_slope_cathode = 0.04
# tafel slope of the anode [V]
tafel_slope_anode = 0.03
# thermo neutral open circuit voltage [V]
v_thermo_neutral = 1.25
# molar concentration of membrane acid groups [mol/m^3]
molar_membrane_acid_group_concentration = 1.2e3
# fitted vapour mass transport coefficient [m/s]
fitted_vapour_vapour_mass_transport_coefficient = 0.62e-5


"""Electric Settings"""
# membrane basic resistance [Ohm/m^2]
membrane_basic_resistance = 0.19
# membrane temperature slope [Ohm/(m^2K)]
membrane_temperature_resistance = 7.e-4
# bipolar plate resistivity [Ohm/m]
bipolar_plate_resistivity = 2.e-6


"""Thermal Settings"""
# thermal conductivity of the bipolar plate through plane  [W/(mK)]
thermal_conductivity_bipolar_plate_z = 1.e2
# thermal conductivity of the bipolar plate in plane  [W/(mK)]
thermal_conductivity_bipolar_plate_x = 1.e2
# thermal conductivity of the gde through plane  [W/(mK)]
thermal_conductivity_gas_diffusion_electrode_z = 1.e0
# thermal conductivity of the gde in plane  [W/(mK)]
thermal_conductivity_gas_diffusion_electrode_x = 1.e0
# thermal conductivity of the membrane through plane  [W/(mK)]
thermal_conductivity_membrane_z = .26e-1
# thermal conductivity of the membrane in plane  [W/(mK)]
thermal_conductivity_membrane_x = .26e-1
# heat capacity of the coolant [J/(kgK)]
heat_capacity_coolant = 2.5e3
# thermal conductivity of the coolant [W/(mK)]
thermal_conductivity_coolant = 0.22
# density of the coolant [kg/m3]
density_coolant = 1052.
# dynamic viscosity of the coolant
dynamic_viscosity_coolant = 56.e-3
# convection coefficient between the stack walls and the environment [W/(Km^2)]
convection_coefficient_stack_environment = 5.

"""Fluid Mechanic Settings"""
# geometrical pressure loss coefficient of the manifold header
manifold_pressure_loss_coefficient = 5.0
# bend pressure loss coefficient of the channel bends
bend_pressure_loss_coefficient = 0.1
# number of channel bends
channel_bends = 48.
# cathode channel gas flow direction
cathode_channel_flow_direction = True
# anode channel gas flow direction
anode_channel_flow_direction = False

"""Humidification"""
# cathode inlet gas relative humidity
inlet_humidity_cathode = 0.0
# anode inlet gas relative humidity
inlet_humidity_anode = 0.0

"""Species settings"""
# oxygen concentration at the cathode inlet
oxygen_inlet_concentration = 0.21
# hydrogen concentration at the anode inlet
hydrogen_inlet_concentration = 0.5

