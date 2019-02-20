ToDo:

**Simon**:
1. Neuen Branch aus dem „develop“ Branch erstellen und verwenden
2. Durch jede Klasse gehen und schauen ob die Member Variablen nötig sind. 
   D.h. die Verwendung von jeder Variable in __init__ die als "self.x" 
   geschrieben wird mit „Find Usages“ in PyCharm checken:
   * Ist diese Variable einzigartig für diese Klasse? 
     Variablen sollten nicht doppelt, dreifach oder vierfach in irgendwelchen
     Objekten gespeichert sein. Bei Integern wie der Elementzahl oder sowas 
     ist es okay, aber keinesfalls bei Arrays? Gehört die Variable nicht zur 
     Klasse, wird aber in der aktuellen Klasse verwendet, kann in der benötigten
     Funktion direkt auf das andere Objekt und seine Variable zugegriffen 
     werden:
   * Wenn sie nur in __init__ verwendet wird:
     okal machen (ohne „self.“)  und direkt an den Ort rücken, wo sie 
     verwendet wird.   
   * Wenn sie nur in 1-2 Funktionen außerhalb von __init__ verwendet wird:
     - Kann sie sinnvoll als input der Funktion verwendet werden? 
     Z.B. wenn es physikalisch eher eine Variable als ein Parameter ist. 
     Denke an Funktion y = a*x^2 + b*x + c 
     (a, b, c sind eher Member Variablen, weil sie nicht so oft variieren, x 
     wird in die Funktion übergeben und y zurückgegeben.)
     - Kann sie aus dem Dictionary gelesen werden, welches oft als input für 
     die __init__ Funktion übergeben wird? Dieses Dictionary muss dann 
     entsprechend als Member Variable gespeichert werden. In der Funktion 
     entweder dann einer lokalen Variable zuordnen oder direkt bei Anwendung 
     aus demDictionary lesen, je nachdem wie oft sie verwendet wird. 
     - Bei Arrays oder Listen sind wahrscheinlich öfter 
     sinnvoll als Member abzuspeichern



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
  
**Species Data**
-   Species data should be stored in a class structure
    including phase as lower level classes.