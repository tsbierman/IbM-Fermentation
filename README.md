# Individual-based Model for Fermentation

*Contributors: Thomas Bierman, Eloi Martinez-Rabert, Chiel van Amstel, Rebeca Gonzalez-Cabaleiro

An individual-based model (IbM), sometimes called an agent-based model (ABM), is a model that considers every individual in the system separately.
In this case, this means that every microorganism will be considered an individual with its own set of metabolic coefficients and kinetic parameters, such as the growth rate (µ). The microorgansisms grow, divide or die independently from each other. The consumption and production rates differ for each organisms as they are based on local concentrations. The substrates and products involved in these reaction diffuse through the system. Approaching the situation with a IbM can capture the heterogenity of the system, regarding both the microorganisms and the concentrations, better than an ODE-approach could. An IbM also considers the interactions between different microorganisms as they are affected by each other's activities.
The individuals are connected to each other through the bulk liquid. The concentrations in the bulk liquid will change depending on the inflow and the outflow of the reactor and the combined activities of all the individuals.

This model is based on the [IbM-framework](https://github.com/Computational-Platform-IbM/IbM) created in another repository. However, that model is based on _Nitrospira_, whereas this model will be based around anaerobic fermentation. The code of this model is mostly written in [Julia](https://julialang.org/), but will also contain parts in Java.
_______________________________

**:warning: To open the links in a new tab: right click on the link + "Open link in new tab". :warning:**

**Methods of IbM can be downloaded [here](https://github.com/Computational-Platform-IbM/IbM/raw/main/Documents/Methods.pdf).**

## :gear: Installations
IbM-Fermentation is build in Julia. Thus, Julia must be installed on your computer. Julia can be downloaded [here](https://julialang.org/downloads/).
<br> Julia is using a command-line Interface. To make working with the code easier, it is also possible to run it with Visual Studio Code (VS Code), which can be downloaded [here](https://code.visualstudio.com/Download). In this case, a Julia extension has to be installed. In VS Code, open "Extensions" by clicking on the blocks on the left or pressing Ctrl+Shift+X. Then simply search for Julia and install the extension (from julialang).
<br> As part of the code is written in Java, a Java Devolopment Kit (JDK) is required. A JDK from Oracle can be downloaded [here](https://www.oracle.com/java/technologies/downloads/).

## :clipboard: Instructions for first time use
1. Download the code as .zip from [here](https://github.com/tsbierman/IbM-Fermentation/releases)
2. Extract the files to a destination (🌟 recommendation: Desktop).
3. Open Julia (or VS Code).
    - For more information about the VS Code User Interface, click [here](https://code.visualstudio.com/docs/getstarted/userinterface).
4. Add the required additional packages. First write `using Pkg`. Then, write `Pkg.add(Package)`, with `Package` being one of the following:
    - XLSX
    - JavaCall
    - Plots
    - StatsPlots
    - InvertedIndices
    - Statistics
    - DSP
    - DifferentialEquations
    - ODE
    - DataStructures
    - FileIO
    - JLD2
    - TickTock
    - ImageFiltering

**:warning: Some of these are quick, others might take a while :warning:**
5. Create a map `results` in the **Code folder<sup>1</sup>**.
<br><sup><sup>1</sup> Code folder: folder with `IbM.jl` file. </sup><br>

## :clipboard: Intruction for model use
1. Open Julia (or VS Code).
2. Go the the **Code folder<sup>2</sup>**.
<br><sup><sup>2</sup> Code folder: folder with `IbM.jl` file. </sup><br>
    → Make sure that `pwd()` yields `~\\IbM-fermentation`, thus the folder that `IbM.jl` is in. <br>
    → Move around with `cd("newFolder")`. <br>
    → Moving a directory upwards can be done with `cd(dirname(pwd()))`. <br>
    → More information about Filesystem commands can be found [here](https://docs.julialang.org/en/v1/base/file/).

3. Create the seed-file<sup>3</sup>.
<br><sup><sup>3</sup> Seed-file: `.jld2` file that stores the variables after initialising. The file is used to execute the code. </sup>
    1. Create a copy of the AF_template (planning\Templates\AF_template.xlsx) and modify the file with all parameters.<br>
    &#09;<sup>Instruction on how to use the Excel can be found in the *Information* sheet.</sup><br>
    2. Place the modified Excel in the "planning" folder
    3. Write `include("inclusion_file.jl")` to *Command Window*.<br>
        - This will initialise a JVM and include all required files for running a simulation. This only needs to be done once per Julia session.<br>
    &#09;<sup>When the model is under development and testing requires functions to be updated after adjustments, disable the initialisation of the JVM in `inclusion_file.jl`.</sup><br>
    4. Write `create_mat("planning\\Excel_Name.xlsx", xxxx);` to *Command Window* <br>(❗where Excel_Name is the name of the Excel and xxxx is the simulation number from 1 to 9999❗).<br>
        - This reads out the excel file and stores the variabels in `sim_xxxx.jld2` in the **Code folder**.

4. Execute IbM code<br>
    1. Check whether the desired seed-file (`sim_xxxx.jld2`) is located in the **Code folder** (folder with `IbM.jl` file).<br>
    2. Check whether a **results** folder is present within Code folder. <br>
    3. Write `IbM(xxxx)` to *Command Window* (❗where xxxx is the simulation number❗).<br>
        - Once the simulation is done, `sim_xxxx.jld2` (i.e. the seed-file) is moved to the corresponding **results** folder.<br>
5. Get Data or Visualisation of Results (see below).
__________________________
## :mag: Additional model information
### Create_mat()
The create_mat function `(lib\\pre_processing\\create_mat.jl)` can function in 2 different ways. The simulation number dictates which one will be executed.<br>
- The normal use case is when the simulation number is between 1 and 9999. `create_mat()` will read out the excel file and save the variables in a `sim_xxxx.jld2`. `IbM(xxxx)` can be called after that. This way, a user can choose to produce several `.jld2` files for several simulations consequetively without each simulation running in between.<br>
- For testing purposes, it is not efficient to save every time `create_mat()` is called. Therefore, when calling the function with a negative simulation number, no `.jld2` file will be created. Instead, `create_mat()` will just return the structs that would have been saved. These can then be used for testing purposes.<br>
__________________________

## Get Data
Results are stored in the **results folder** (`results\xxxx` where `xxxx` is the simulation number) in the `results1D.jld2` or `results2D.jld2` file. In all cases, first navigate to the **Code folder** and write `include("inclusion_file.jl")`.

If the data is required in the Julia environment, it can be accessed using `variable_name = load(results_file_name, "variable_name")`, where `results_file_name` is the path to the result file and `variable_name` is the variable that needs to be extracted. Multiple variables can be extracted at the same time. Which variables are saved can be seen in the saving scripts (in `lib\\post_processing`). Modifying what is saved can be done there as well.

If the data is required in an Excel file, the package [XLSX](https://github.com/felipenoris/XLSX.jl) is used. First the data has to be loaded into the environment, as described above. Then XLSX is used to write the data to an Excel file. [This](https://felipenoris.github.io/XLSX.jl/stable/tutorial/) XLSX tutorial shows how this package is used. Additionally, the package is used in the `extract_data.jl` script, which can be used as reference.

__________________________

## Visualize data
To start, navigate to the **Code folder** and write `include("inclusion_file.jl")`. This also makes all visualisation scipts available. Information about the input and output of these scripts can be found in the scripts themselves (in the folder `lib\\post_processing`). 
