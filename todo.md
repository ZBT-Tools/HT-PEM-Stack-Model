ToDo:

**Manifold:**
-   The viscosity of the gas-mixtures in the outlet manifold should be 
    calculated using Herning and Zipperer. At the moment the viscosity of the
    gas-mixture along the outlet header is equal to each local fuel cell 
    outlet viscosity. To calculate the gas mixture viscosity over the outlet 
    manifold, it is necessary to calculate the species concentrations in the 
    outlet manifold.

-   The thermal coupling between the manifold-channels and the bipolar-plate 
    should be added. At the moment the manifold-channels are idealized as 
    adiabatic. Therefore it is necessary to modify the thermodynamic system.
 
-   The effects of fluid water in the outlet manifolds should be added. 
    Direct effects of fluid water in the outlet manifolds are neglected at 
    the moment.
 
**Electrical Coupling**
-   The effects of under stoichiometry are not completely implemented yet.
 
**Thermal System**
-   The coolant channel and species channels are implemented explicitly. To 
    decrease the number of necessary iterations to solve the thermal system, 
    the coolant and species channel could be implemented implicitly.
-   Another option would be to decouple the implicit calculation of thermal 
    conduction between all fuel cells. Instead of solving one large thermal 
    matrix, multiple smaller matrices would be solved individually and
    additional iterations would be needed. If the number of iterations for 
    the decoupled system stays rather small, this could be beneficial for 
    large stacks. The number of iterations should be small assuming that the 
    temperature differences between adjacent bipolar plates are small. 
  
**Polarization Curve**
-   Some automatization algorithm to fit given polarization curves
    using defined fitting-parameter should be implemented.
