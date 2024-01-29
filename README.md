# IbM-Fermentation

*Contributors: Thomas Bierman

An individual-based model (IbM), sometimes called an agent-based model (ABM), is a model that considers every individual in the system separately.
In this case, this means that every micro-organism will be considered an individual with its own set of metabolic coefficients and kinetic parameters, such as the growth rate µ. The micro-organsism thus grow, divide or die independently from each other. The consumption and production rates differ for each organisms as they are based on local concentrations. The substrates and products involved in these reaction diffuse through the system. Approaching the system with a IbM can capture the heterogenity of the system, regarding both the micro-organisms and the concentrations, better than an ODE-approach could. An IbM also considers the interactions between different micro-organisms as they are affected by each other's activities.
The individuals are connected to each other through the bulk liquid. The concentrations in the bulk liquid will change depending on the inflow of the reactor and the combined activities of all the individuals.

This model is based on the [IbM-framework](https://github.com/Computational-Platform-IbM/IbM) created in another repository. However, that model is based on _Nitrospira_, whereas this model will be based around
anaerobic fermentation of glucose. The code of this model is mostly written in [Julia](https://julialang.org/), but will also contain parts in Java.
_______________________________

**:warning: To open the links in a new tab: right click on the link + "Open link in new tab". :warning:**

**Methods of IbM can be downloaded [TODO]().**

## :gear: Julia installation
IbM-Fermentation is build up in Julia. Thus, Julia must be installend on your computer.
<br> Julia can be downloaded [here](https://julialang.org/downloads/).
<br> Julia is using a command-line Interface. To make working with the code easier, it is also possible to run it with Visual Studio Code (VS Code), which can be downloaded [here](https://code.visualstudio.com/Download). In this case, a Julia extension has to be installed. In VS Code, open "Extensions" by clicking on the blocks on the left or pressing Ctrl+Shift+X. Then simply search for Julia and install the extension (from julialang).

## :clipboard: Instructions for model use
1. Download the code as .zip. Last version **TODO** [Download code]().
2. Extract the files to a destination (🌟 recommendation: Desktop).
3. Open Julia (or VS Code).
    - For more information about the VS Code User Interface, click [here](https://code.visualstudio.com/docs/getstarted/userinterface).
4. Go the the **Code folder<sup>1</sup>**.
<br><sup><sup>1</sup> Code folder: folder with `IbM.jl` file. </sup><br>
    → Make sure that `pwd()` yields `~\\IbM-fermentation`, thus the folder that `IbM.jl` is in. <br>
    → Move around with `cd("newFolder")`. <br>
    → Moving a directory upwards can be done with `cd(dirname(pwd()))`. <br>
    → More information about Filesystem commands can be found [here](https://docs.julialang.org/en/v1/base/file/).

5. Create the seed-file<sup>2</sup>.
<br><sup><sup>2</sup> Seed-file: `.jld2` file that stores the variables after initialising. The file is used to execute the code. </sup>
    1. Modify the main Excel (lib\planning\Excels\main.xlsx) with all parameters.<br>
    &#09;<sup>Instruction on how to use main.xlsx can be found in *Information* sheet.</sup><br>
    2. Write `include("lib\\pre_processing\\initialiseJVM.jl")` to *Command Window*.<br>
        - This initiates the Java Virtual Machine (JVM). This only needs to be done once per Julia session.<br>
    3. Write `include("inclusion_file.jl")` to *Command Window*.<br>
        - This will import all required modules and include all required files for running a simulation. This only needs to be done once per Julia session.<br>
    4. Write `create_mat("planning\\main.xlsx", xxxx)` to *Command Window* <br>(❗where xxxx is the simulation number from 1 to 9999❗).<br>
        - This reads out the excel file and stores the variabels in `sim_xxxx.jld2` in the Code folder.

6. Execute IbM code<br>
    1. Check whether the desired seed-file (`sim_xxxx.jld2`) is located in the Code folder (folder with `IbM.jl` file).<br>
    2. Write  `IbM(xxxx)` to *Command Window* (❗where xxxx is the simulation number❗).<br>
        - Once the simulation is done, `sim_xxxx.jld2` (i.e. the seed-file) is moved to the corresponding **results** folder.<br>
7. Get Data or Visualisation of Results (TODO: add instructions later).
__________________________
## :mag: Additional model information
### Create_mat()
The create_mat function `(lib\\pre_processing\\create_mat.jl)` can function in 2 different ways. The simulation number dictates which one will be executed.<br>
- The normal use case is when the simulation number is between 1 and 9999. `create_mat()` will read out the excel file and save the variables in a `sim_xxxx.jld2`. The user will have to call the `IbM(xxxx)` themselves. This way, a user can choose to produce several `.jld2` files for several simulations consequetively without each simulation running inbetween.<br>
- For testing purposes, it is not efficient to save every time `create_mat()` is called. Therefore, when calling the function with a negative simulation number, no `.jld2` file will be created. Instead, `create_mat()` will just return the structs that would have been saved. These can then be used for testing properties.<br>
