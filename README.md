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
1. Download the code as .zip. Last version **TODO**. [Download code]()
2. Extract the files to a destination (🌟 recommendation: Desktop)
3. Open Julia (or VS Code).
    <br>&emsp;- For more information about the VS Code User Interface, click [here](https://code.visualstudio.com/docs/getstarted/userinterface)
4. Go the the **Code folder<sup>1</sup>**
    <br><sup><sup>1</sup> Code folder: folder with `IbM.jl` file. </sup>
    <br>→ Make sure that `pwd()` yields `~\\IbM-fermentation`, thus the folder that `IbM.jl` is in
    <br>&emsp;→ Move around with `cd("newFolder")`
    <br>&emsp;→ Moving a directory upwards can be done with `cd(dirname(pwd()))`
    <br>&emsp;→ More information about Filesystem commands can be found [here](https://docs.julialang.org/en/v1/base/file/).
5. Create the seed-file<sup>2</sup>
    <br><sup><sup>2</sup> Seed-file: `.jld2` file that stores the variables after initialising. The file is used to execute the code </sup>
    <br>&emsp;1. Modify the main Excel (lib\planning\Excels\main.xlsx) with all parameters.
        <br><sup>Instruction on how to use main.xlsx can be found in *Information* sheet.</sup>
    <br>&emsp;2. Write `include("lib\\pre_processing\\initialiseJVM.jl")` to *Command Window*
        <br>- This initiates the Java Virtual Machine (JVM). This only needs to be done once per Julia session.
    <br>&emsp;3. Write `include("lib\\pre_processing\\create_mat.jl");` to *Command Window*
        <br>&emsp;&emsp;- This will run a complete workflow with `start_up.xlsx`. This excel will be loaded out and a short simulation will be run with the saved data. This could take a while (TODO:time_indication), because Julia will have to compile everything. However, this compiling will speed up the actual simulation. For optimal simulation time, this should be done once every Julia session.
        &emsp;<br>- This will create a map with the name `0000` in **results**
    <br>&emsp;4. Write `create_mat("planning\\main.xlsx", xxxx)` to *Command Window* (❗ where xxxx is the simulation number from 1 to 9999)
        <br>&emsp;&emsp;- This reads out the excel file and stores the variabels in `sim_xxxx.jld2` in the Code folder.
6. Execute IbM code:
    <br>&emsp;1. Check whether the desired seed-file (`sim_xxxx.jld2`) is located in the Code folder (folder with `IbM.jl` file). 
    <br>&emsp;2. Write `include("IbM.jl")` to *Command Window*
    <br>&emsp;3. Write  `IbM(xxxx)` to *Command Window* (❗ where xxxx is the simulation number)
    <br>&emsp;&emsp;- Once the simulation is done, `sim_xxxx.jld2` (i.e. the seed-file) is moved to the corresponding **results** folder.
7. Get Data or Visualisation of Results (TODO: add instructions later)


