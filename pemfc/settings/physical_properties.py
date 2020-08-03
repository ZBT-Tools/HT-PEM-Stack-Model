""" Physical properties
(i.e. material properties, electrochemical settings, etc.)"""

"""Electrochemistry Settings"""
# volumetric exchange current density of the cathode[A/m^3]
exchange_current_density_cathode = 8.0e4
# volumetric exchange current density of the anode [A/m^3]
exchange_current_density_anode = 5.0e8
# oxygen catalyst layer diffusion coefficient [m^2/s]
oxygen_catalyst_layer_diffusion_coefficient = 5e-8
# hydrogen catalyst layer diffusion coefficient [m^2/s]
hydrogen_catalyst_layer_diffusion_coefficient = 1e-7
# oxygen gas diffusion layer diffusion coefficient [[m^2/s]
oxygen_gas_diffusion_layer_diffusion_coefficient = 6e-6  # 2.95e-6
# hydrogen gas diffusion layer diffusion coefficient [m^2/s]
hydrogen_diffusion_layer_diffusion_coefficient = 10.0e-6
# catalyst layer proton conductivity of the cathode [Ohm^-1/m]
catalyst_layer_proton_conductivity_cathode = 1.2
# catalyst layer proton conductivity of the anode [Ohm^-1/m]
catalyst_layer_proton_conductivity_anode = 1.2
# tafel slope of the cathode [V]
tafel_slope_cathode = 0.035
# tafel slope of the anode [V]
tafel_slope_anode = 0.035
# thermo neutral open circuit voltage [V]
v_thermo_neutral = 1.25

"""Membrane Resistance"""
# membrane type ('Constant', 'Springer', or 'Kvesic') (Springer-Model
# implementation must still be checked properly)
membrane_type = 'Kvesic'
# membrane basic resistance [Ohm-m2]
membrane_basic_resistance = 4.3e-5
# membrane temperature slope [Ohm-m2/K)]
membrane_temperature_coefficient = 7e-8
# membrane basic conductivity [S/m] (for 'Constant'-Model)
membrane_basic_conductivity = 5.0
# bipolar plate conductivity [S/m]
electrical_conductivity_bipolar_plate = 20.0
# gas diffusion electrode material conductivity [S/m]
electrical_conductivity_gde = 500.0

"""Membrane Water Transport Properties"""
# molar concentration of membrane acid groups [mol/m^3]
molar_membrane_acid_group_concentration = 1.2e3
# fitted vapour mass transport coefficient [m/s]
vapour_mass_transport_coefficient = 0.62e-5

"""Thermal Settings"""
# thermal conductivity of the bipolar plate through-plane  [W/(mK)]
thermal_conductivity_bipolar_plate_z = 1.e2
# thermal conductivity of the bipolar plate in-plane  [W/(mK)]
thermal_conductivity_bipolar_plate_x = 1.e2
# thermal conductivity of the gde through-plane  [W/(mK)]
thermal_conductivity_gas_diffusion_electrode_z = 1.e0
# thermal conductivity of the gde in-plane  [W/(mK)]
thermal_conductivity_gas_diffusion_electrode_x = 1.e0
# thermal conductivity of the membrane through-plane  [W/(mK)]
thermal_conductivity_membrane_z = .26e0
# thermal conductivity of the membrane in-plane  [W/(mK)]
thermal_conductivity_membrane_x = .26e0

# constant coolant properties
coolant_name = 'FRAGOLTHERM S-15-A at 60°C'
# specific heat of the coolant [J/(kg-K)]
heat_capacity_coolant = 2010.0
# thermal conductivity of the coolant [W/(mK)]
thermal_conductivity_coolant = 0.159
# density of the coolant [kg/m3]
density_coolant = 981.0
# dynamic viscosity of the coolant [Pa-s]
dynamic_viscosity_coolant = 9.81e-3



