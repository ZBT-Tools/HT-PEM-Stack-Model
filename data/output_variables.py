cell_var = [["cathode.mol_flow[0]", "Molar_Flow_", "$mol/s$"],
            ["cathode.mol_flow[1]", "Molar_Cathode_Water_Flow", "$mol/s$"],
            ["cathode.mol_flow[2]", "Molar__Cathode_Nitrogen_Flow", "$mol/s$"],
            ["anode.mol_flow[0]", "Molar_Hydrogen_Flow", "$mol/s$"],
            ["anode.mol_flow[1]", "Molar_Anode_Water_Flow", "$mol/s$"],
            ["anode.mol_flow[2]", "Molar_Anode_Nitrogen_Flow", "$mol/s$"],
            ["cathode.gas_con[0]", "Molar_Oxygen_Concentration", "$mol/m³$"],
            ["cathode.gas_con[1]", "Molar_Cathode_Water_Gas_Phase_Concentration", "$mol/m³$"],
            ["cathode.gas_con[2]", "Molar__Cathode_Nitrogen_Concentration", "$mol/m³$"],
            ["anode.gas_con[0]", "Molar_Hydrogen_Concentration", "$mol/m³$"],
            ["anode.gas_con[1]", "Molar_Anode_Water_Gas_Phase_Concentration", "$mol/m³$"],
            ["anode.gas_con[2]", "Molar_Anode_Nitrogen_Concentration", "$mol/m³$"],
            ["cathode.mass_f[0]", "Oxygen_Mass_Fraction", "%"],
            ["cathode.mass_f[1]", "Cathode_Water_Mass_Fraction", "%"],
            ["cathode.mass_f[2]", "Cathode_Nitrogen_Mass_Fraction", "%"],
            ["anode.mass_f[0]", "Hydrogen_Mass_Fraction", "%"],
            ["anode.mass_f[1]", "Anode_Water_Mass_Fraction", "%"],
            ["anode.mass_f[2]", "Anode_Nitrogen_Mass_Fraction", "%"],
            ["cathode.mol_f[0]", "Oxygen_Mol_Fraction", "%"],
            ["cathode.mol_f[1]", "Cathode_Water_Mol_Fraction", "%"],
            ["cathode.mol_f[2]", "Cathode_Nitrogen_Mol_Fraction", "%"],
            ["anode.mol_f[0]", "Hydrogen_Mol_Fraction", "%"],
            ["anode.mol_f[1]", "Anode_Water_Mol_Fraction", "%"],
            ["anode.mol_f[2]", "Anode_Nitrogen_Mol_Fraction", "%"],
            ["cathode.cp[0]", "Oxygen_Heat_Capacity", "$J/kgK$"],
            ["cathode.cp[1]", "Cathode_Water_Heat_Capacity", "$J/kgK$"],
            ["cathode.cp[2]", "Cathode_Nitrogen_Heat_Capacity", "$J/kgK$"],
            ["anode.cp[0]", "Hydrogen_Heat_Capacity", "$J/kgK$"],
            ["anode.cp[1]", "Anode_Water_Heat_Capacity", "$J/kgK$"],
            ["anode.cp[2]", "Anode_Nitrogen_Heat_Capacity", "$J/kgK$"],
            ["cathode.lambdas[0]", "Oxygen_Heat_Conductivity", "$W/mK$"],
            ["cathode.lambdas[1]", "Cathode_Water_Heat_Conductivity", "$W/mK$"],
            ["cathode.lambdas[2]", "Cathode_Nitrogen_Heat_Heat_Conductivity", "$W/mK$"],
            ["anode.lambdas[0]", "Hydrogen_Heat_Heat_Conductivity", "$W/mK$"],
            ["anode.lambdas[1]", "Anode_Water_Heat_Conductivity", "$W/mK$"],
            ["anode.lambdas[2]", "Anode_Nitrogen_Heat_Conductivity", "$W/mK$"],
            ["cathode.visc[0]", "Oxygen_Viscosity", "$kg/ms$"],
            ["cathode.visc[1]", "Cathode_Water_Viscosity", "$kg/ms$"],
            ["cathode.visc[2]", "Cathode_Nitrogen_Viscosity", "$kg/ms$"],
            ["anode.visc[0]", "Hydrogen_Viscosity", "$kg/ms$"],
            ["anode.visc[1]", "Anode_Water_Viscosity", "$kg/ms$"],
            ["anode.visc[2]", "Anode_Nitrogen_Viscosity", "$kg/ms$"],
            ["cathode.r_gas", "Cathode_Gas_Constant", "$J/kgK$"],
            ["anode.r_gas", "Cathode_Gas_Constant", "$J/kgK$"],
            ["cathode.cp_gas", "Cathode_Heat_Capacity", "$J/kgK$"],
            ["anode.cp_gas", "Anode_Heat_Capacity", "$J/kgK$"],
            ["cathode.lambda_gas", "Cathode_Heat_Conductivity", "$W/mK$"],
            ["anode.lambda_gas", "Anode_Heat_Conductivity", "$W/mK$"],
            ["cathode.rho_gas", "Cathode_Gas_Density", "$kg/m³$"],
            ["anode.rho_gas", "Anode_Gas_Density", "$kg/m³$"],
            ["cathode.u", "Cathode_Flow_Velocity", "$m/s$"],
            ["anode.u", "Anode_Flow_Velocity", "$m/s$"],
            ["cathode.p", "Cathode_Channel_Pressure", "$Pa$"],
            ["anode.p", "Anode_Channel_Pressure", "$Pa$"],
            ["cathode.cp_fluid", "Cathode_Fluid_Water_Heat_Capacity", "$J/kgK$"],
            ["anode.cp_fluid", "Anode_Fluid_Water_Heat_Capacity", "$J/kgK$"],
            ["cathode.m_flow_fluid", "Cathode_Mass_Flow", "$kg/s$"],
            ["anode.m_flow_fluid", "Anode_Mass_Flow", "$kg/s$"],
            ["cathode.temp_fluid", "Cathode_Fluid_Temperature", "$K$"],
            ["cathode.act_loss", "Cathode_Activation_Overpotential", "$V$"],
            ["anode.act_loss", "Anode_Activation_Overpotential", "$V$"],
            ["cathode.cl_diff_loss", "Cathode_Catalyst_Layer_Diffusion_Overpotential", "$V$"],
            ["anode.cl_diff_loss", "Anode_Catalyst_Layer_Diffusion_Overpotential", "$V$"],
            ["cathode.gdl_diff_loss", "Cathode_Gas_Diffusion_Layer_Overpotential", "$V$"],
            ["anode.gdl_diff_loss", "Anode_Gas_Diffusion_Layer_Overpotential", "$V$"],
            ["mem_loss", "Membran_Overpotential", "$V$"],
            ["v_loss", "Voltages_Losses", "$V$"],
            ["v", "Cell_Voltage", "$V$"]]

stack_var = [["temp_cpl_stack.temp_cool", "Coolant_Fluid_Temperature", "K"],
            ["i_ca", "Current_Density", "A/m²"],
            ["manifold[0].cell_stoi", "Cathode Stoichiometry", ""],
            ["manifold[1].cell_stoi", "Anode Stoichiometry", ""]]






