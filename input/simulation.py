""" This file contains the simulation settings"""

"""Simulation Settings"""
# partition of the flow channel in elements along the x-axis
elements = 10
# convergence criteria of the simulation
convergence_criteria = 1.e-6
# maximal number of iterations
maximal_number_iteration = 50
# output csv data
save_csv_data = True
# output plots
save_plot_data = True
# calculate the PEMFC stack temperatures
calc_temperature = True
# calculate the current density distribution
calc_current_density = True
# calculate the flow distribution
calc_flow_distribution = True
# calculate the activation voltage losses
calc_activation_loss = True
# calculate the membrane voltage losses
calc_membrane_loss = True
# calculate gdl diffusion voltage losses
calc_gdl_loss = True
# calculate catalyst layer diffusion voltage losses
calc_cl_loss = True
# show voltage losses in the voltage-current-density-graph
show_voltage_loss = False

