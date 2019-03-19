To-do-list:
-
The development of the PEMFC-Stack-Model is ist in progress.

The main issues can be found in the manifold,
 thermic and electrical couplings are:

**Manifold:**
- The viscosity of the gas-mixtures in the outlet manifold should be calculated
  using Herning and Zipperer.  
  At the moment the viscosity of the gas-mixture
  of the in the t-junction ending cell is taken
  to be the absolute viscosity at this point.
  To calculated the gas mixture viscosity over the outlet manifold
  it is necessary to calculate all species flows over the outlet manifold.

- The thermic coupling between the manifold-channels
  and the bipolar-plate should be added.
  At the moment the manifold-channels are idealized as adiabatic.
  Here fore it is necessary to modify the thermodynamic system.
 
- The effects of fluid water in the outlet manifolds should be added.
  Direct effects of fluid water in the outlet manifolds are neglected
  at the moment.
 
 **Electrical Coupling**
 - The effects of under stoichiometry are not completely implemented yet.
 
 **Thermic System**
- The coolant channel and species channels are implemented explicit.
  To decrease the number of necessary iterations to solve the
  thermal system the coolant and species channel should be implemented implicit.
- To increase the calculation speed for larger stacks with high cell number and
  high grid density, decreasing the temperature system level down
  to the cell level is necessary. This should also ease the understanding of the
  matrix position of connecting heat conductivities.

  
**Polarization Curve**
- Some algorithm to fit given polarization curves
  using defined fitting-parameter should be implemented.
  
**Species Data**
- Species data should be stored in an class structure
 including phase as lower level classes.
 
**Classes and Methods**
- Methods might be set as stand alone classes: For example the Class VoltageLoss
  including the fix calculations, inside the init definition and the dynamic
  calculations as methods.
  This would either clear the init space of the higher class levels like
  HalfCell, but also create a more suitable workspace.
  
**Channel**
- It should be considered to treat the coolant channel like the gas channels,
  that means:
  - examine the coolant distribution over the stack
  - examine the pressure drop in the coolant channel
  - examine the species dependence of the coolant on temperature and pressure
    here for some specific data is needed
    
  An option here is to rebuild the class channel on an universal fundament,
  including flow mechanics.

    
   

 
       