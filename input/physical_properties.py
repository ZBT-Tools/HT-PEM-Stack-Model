""" This file contains the physical settings"""

"""Electrochemistry Settings"""
# volumetric exchange current density of the cathode[A/m^3]
exchange_current_density_cathode = 0.75e6
# volumetric exchange current density of the anode [A/m^3]
exchange_current_density_anode = 0.817e9
# oxygen catalyst layer diffusion coefficient [m^2/s]
oxygen_catalyst_layer_diffusion_coefficient = 0.20e-7
# hydrogen catalyst layer diffusion coefficient [m^2/s]
hydrogen_catalyst_layer_diffusion_coefficient = 5.e-8
# oxygen gas diffusion layer diffusion coefficient [[m^2/s]
oxygen_gas_diffusion_layer_diffusion_coefficient = 6.75e-6
# hydrogen gas diffusion layer diffusion coefficient [m^2/s]
hydrogen_diffusion_layer_diffusion_coefficient = 9.52e-6
# catalyst layer proton conductivity of the cathode [Ohm^-1/m]
catalyst_layer_proton_conductivity_cathode = 3.e0
# catalyst layer proton conductivity of the anode [Ohm^-1/m]
catalyst_layer_proton_conductivity_anode = 3.e0
# tafel slope of the cathode [V]
tafel_slope_cathode = 0.03
# tafel slope of the anode [V]
tafel_slope_anode = 0.03
# thermo neutral open circuit voltage [V]
v_thermo_neutral = 1.25
# molar concentration of membrane acid groups [mol/m^3]
molar_membrane_acid_group_concentration = 1.2e3
# fitted vapour mass transport coefficient [m/s]
vapour_mass_transport_coefficient = 0.62e-5

"""Electric Settings"""
# membrane type
membrane_type = 'Springer'
# membrane basic resistance [Ohm/m^2]
membrane_basic_resistance = 0.19
# membrane basic conductivity [S/m]
membrane_basic_conductivity = 10.0
# membrane temperature slope [Ohm/(m^2K)]
membrane_temperature_coefficient = 7.e-4
# bipolar plate conductivity [S/m]
electrical_conductivity_bipolar_plate = 30.0
# gas diffusion electrode material conductivity [S/m]
electrical_conductivity_gde = 5000.0

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
thermal_conductivity_membrane_z = .26e0
# thermal conductivity of the membrane in plane  [W/(mK)]
thermal_conductivity_membrane_x = .26e0
heat_capacity_coolant = 2.5e3
# thermal conductivity of the coolant [W/(mK)]
thermal_conductivity_coolant = 0.22
# density of the coolant [kg/m3]
density_coolant = 1052.
# dynamic viscosity of the coolant
dynamic_viscosity_coolant = 56.e-3



