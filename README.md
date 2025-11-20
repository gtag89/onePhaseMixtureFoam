# onePhaseMixtureFoam

About onePhaseMixtureFoam
-------------------------

onePhaseMixtureFoam is a free, open-source CFD solver based on OpenFOAM capable of simulating multiphase flow instabilities in unsaturated porous media. The solver employs a hyperbolic momentum formulation derived from mixture theory to capture non-equilibrium phenomena such as saturation overshoot and gravity fingering that cannot be reproduced by conventional Richards equation-based approaches.

The code is developed by Dr. Giorgia Tagliavini (https://github.com/gtagliavini - Laboratory for Computational Engineering, Empa, Swiss Federal Laboratories for Materials Science and Technology, Dübendorf, Switzerland).

The implementation is based on mixture theory with a simplified single-phase formulation that neglects the density and viscosity of air. A second viscosity term transforms the momentum equation from parabolic to hyperbolic, enabling finite-speed traveling waves and non-equilibrium front dynamics. The mathematical formulation was conceptualized by Dr. Ivan Lunati (Laboratory for Computational Engineering, Empa, Swiss Federal Laboratories for Materials Science and Technology, Dübendorf, Switzerland).

Code capabilities
-----------------
* Modified PIMPLE algorithm for saturation-velocity coupling instead of traditional pressure-velocity coupling
* Van Genuchten-Mualem capillary pressure relationships with hysteresis framework
* Second viscosity implementation with semi-implicit correction for numerical robustness  
* Captures saturation overshoot across five orders of magnitude in infiltration rates
* Reproduces gravity fingering instabilities in 2D and 3D heterogeneous media
* Excellent parallel scaling performance 

Compatibility
-------------
The code is prepared for compilation with OpenFOAM v2212 and later versions (https://openfoam.org)

Compilation
-----------
Note: the scripts must be run from terminals with sourced OpenFOAM v2212+

```bash
# Navigate to solver directory
cd onePhaseMixtureFoam/solver

# Compile the solver
wmake

# Verify installation
which onePhaseMixtureFoam
```

Tutorials
---------
* **1D_overshoot**: Demonstrates saturation overshoot formation during vertical infiltration, validated against DiCarlo (2004) experiments

![Fig_09](https://github.com/user-attachments/assets/4a18bde5-3c2f-42ca-b6df-88d8093c9f04)

* **2D_fingers**: Shows gravity fingering instabilities in heterogeneous porous media with random permeability distribution

  ![Fig_10](https://github.com/user-attachments/assets/c89b8a57-5452-47ee-a613-234b76083acd)

Each tutorial directory contains:
- Complete OpenFOAM case setup (0/, constant/, system/)
- ./run script with case description 

For users
---------
Please note that this solver is developed by a small research team. The onePhaseMixtureFoam solver is distributed, similarly to OpenFOAM itself, in the hope that it will be useful but without any warranty. We welcome feedback and bug reports through GitHub Issues.

For questions about specific applications or solver capabilities, please use GitHub Discussions. If you need to contact the authors directly regarding onePhaseMixtureFoam, please email: giorgia.tagliavini89@gmail.com .

Cite this work as
-----------------
* Tagliavini, G.; Lunati, I.: onePhaseMixtureFoam: A hyperbolic momentum solver for saturation overshoot and gravity fingering in unsaturated porous media. [In preparation for: Advances in Water Resources].

License
-------
onePhaseMixtureFoam is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.

****************************************************************************
