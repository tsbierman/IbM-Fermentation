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
2. Extract the files to a destination (:bulb: recommendation: Desktop)
3. Open Julia (or VS Code).
    - For more information about the VS Code User Interface, click [here](https://code.visualstudio.com/docs/getstarted/userinterface)
4. Go the the **Code folder<sup>1</sup>**
    <br><sup><sup>1</sup> Code folder: folder with `IbM.jl` file. </sup>
    → Make sure that `pwd()` yields `~\\IbM-fermentation`, thus the folder that `IbM.jl` is in
    → Move around with `cd("newFolder")`
    → Moving a directory back up can be done with `cd(dirname(pwd()))`
    → More information about Filesystem commands can be found [here](https://docs.julialang.org/en/v1/base/file/)
5. Create the seed-file<sup>2</sup>
    <br><sup><sup>2</sup> Seed-file: `.jld2` file that stores the variables after initialising. The file is used to execute the code </sup>
    → Write `include("lib\\pre_processing\\initialiseJVM.jl")`
        - This initiates the Java Virtual Machine (JVM)
        - This only needs to be done once per session
    → Modify the main Excel (lib\planning\Excels\main.xlsx) with all parameters.
      <br><sup>Instruction on how to use main.xlsx can be found in *Information* sheet.



THOUGHTS:
INCLUDING CREATE_MAT:
include("lib\\pre_processing\\create_mat.jl");
this should check whether this is run like this, if so, run a 0000 file that will automatically be sent into the IbM, this way, everything has been loaded and actually simulation should take shorter.

RUNNING CREATE_MAT example:
create_mat("planning\\test_file.xlsx");

SAVING NEEDS TO BE DONE WITHIN THE CREATE_MAT


THEN RUN THE IbM code 
CHECK WHETHER A FILE CALLED XXXX.jlD2 IS IN SAME DIRECTORY AS IbM.jl
    If not move it there
IbM(XXXX)
AFTERWARDS .jld2 WILL BE MOVED TO THE RESULTS FOLDER