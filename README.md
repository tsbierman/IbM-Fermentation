# IbM-Fermentation
---
An individual-based model (IbM) sometimes called an agent-based model (ABM), is a model that consideres every individual in the system separately.
In this case, this means that every micro-organism will be considered an individual with its own set of metabolic reaction and kinetic parameters.
Therefore, when running such a model, the growth, decay, uptake and production rates will be calculated for every individual separately.
The individuals are connected to each other by the bulk liquid. The concentrations in the bulk liquid will change depending on the 
inflow of the reactor and the combined activities of all the individuals.

This model is based on the ([IbM](https://github.com/Computational-Platform-IbM/IbM)) created in another repository. However, that model is based on _Nitrospira_, whereas this model will be based around
anaerobic fermentation of glucose. The code of this model is mostly written in ([Julia](https://julialang.org/)).
