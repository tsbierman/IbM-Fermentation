# IbM-Fermentation

*Contributors: Thomas Bierman

An individual-based model (IbM), sometimes called an agent-based model (ABM), is a model that considers every individual in the system separately.
In this case, this means that every micro-organism will be considered an individual with its own set of metabolic coefficients and kinetic parameters, such as the growth rate µ. The micro-organsism thus grow, divide or die independently from each other. The consumption and production rates differ for each organisms as they are based on local concentrations. The substrates and products involved in these reaction diffuse through the system. Approaching the system with a IbM can capture the heterogenity of the system, regarding both the micro-organisms and the concentrations, better than an ODE-approach could. An IbM also considers the interactions between different micro-organisms as they are affected by each other's activities.
The individuals are connected to each other through the bulk liquid. The concentrations in the bulk liquid will change depending on the inflow of the reactor and the combined activities of all the individuals.

This model is based on the [IbM-framework](https://github.com/Computational-Platform-IbM/IbM) created in another repository. However, that model is based on _Nitrospira_, whereas this model will be based around
anaerobic fermentation of glucose. The code of this model is mostly written in [Julia](https://julialang.org/), but will also contain parts in Java.
_______________________________

## Julia installation
IbM-Fermentation is build up in Julia. Thus, Julia must be installend on your computer.
<br> Julia can be downloaded [here](https://julialang.org/downloads/).
<br> Julia is written from a terminal. To make working with the code easier, it is also possible to run it with Visual Studio Code (VSC), which can be downloaded [here](https://code.visualstudio.com/Download). In this case, a Julia extension has to be installed. In VSC, open extension by clicking on the blocks on the left or pressing Ctrl+Shift+X. Then simply search for Julia and install the extension (from julialang).

## Usage instructions
