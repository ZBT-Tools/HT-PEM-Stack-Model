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
 -  The coolant channel and species channels are implemented explicit.
  To decrease the number of necessary iterations to solve the
  thermal system the coolant and species channel should be implemented implicit.
  - Decreasing the temperature-matrix to cell level might increase the performance.

  
**Polarization Curve**
- Some algorithm to fit given polarization curves
  using defined fitting-parameter should be implemented.
  
**Species Data**
- Species data should be stored in an class structure
 including phase as lower level classes.
