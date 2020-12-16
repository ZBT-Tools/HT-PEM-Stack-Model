# PEMFC-Stack-Model
A reduced dimensional numerical model to simulate the performance of PEM fuel cell stacks developed in Python 3.6 utilizing the numerical libraries NumPy and SciPy.

Features:
- physical stack domain is discretized into two dimensions:
    - through each cell in the direction of the electrical current (current-direction)
    - along the flow direction of each channel (flow-direction)

- calculation of the reactant flow distribution into the cells 
  based on the geometry of headers and channels
  
- local current distribution along the flow- and current-direction due to:
    - reactant transport within the channels and the porous media
    - temperature distribution
    - reaction kinetics and voltage losses according to Kulikovsky (2013)
     
- temperature distribution along the flow- and current-direction with a discretization in the current-direction (through plane) in five nodes at the interfaces of:
    - anodic and cathodic bipolar plates (BPP-BPP)
    - anodic bipolar plate and gas diffusion electrode (BPP-GDE, Ano)
    - anodic gas diffusion electrode and membrane (GDE-Mem, Ano)
    - cathodic gas diffusion electrode and membrane (GDE-Mem, Cat)        
    - cathodic bipolar plate and gas diffusion electrode (BPP-GDE, Cat)

# Minimum requirements:
- NumPy 1.14.3
- SciPy 1.1.0
- Matplotlib 2.2.2

# Usage
Download the repository and execute the pemfc/main_app.py (pure cli program) or 
the pemfc/gui_app.py (GUI program) file with your Python interpreter. Input
parameters can be adapted in the corresponding files in the "settings" folder. 
At the end of a simulation run, a folder called "output" will be created, 
which contains the results in various data files and plots. Settings can be
 adjusted in the GUI or in the specific files under pemfc/settings.

# References:
Stack discretization, temperature coupling, reactant transport and membrane properties according to:  
*Chang, Paul, Gwang-Soo Kim, Keith Promislow, und Brian Wetton. „Reduced Dimensional Computational Models of Polymer Electrolyte Membrane Fuel Cell Stacks“. Journal of Computational Physics 223, Nr. 2 (Mai 2007): 797–821. https://doi.org/10.1016/j.jcp.2006.10.011.*

Manifold model and flow distribution calculation according to:  
*Koh, Joon-Ho, Hai-Kyung Seo, Choong Gon Lee, Young-Sung Yoo, und Hee Chun Lim. „Pressure and flow distribution in internal gas manifolds of a fuel-cell stack“. Journal of Power Sources 115, Nr. 1 (2003): 54–65.*

Electrochemical reaction kinetics and transport losses according to:  
*Kulikovsky, A. A. „A Physically-Based Analytical Polarization Curve of a PEM Fuel Cell“. Journal of the Electrochemical Society 161, Nr. 3 (28. Dezember 2013): F263–70. https://doi.org/10.1149/2.028403jes.*


