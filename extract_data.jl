function extract_data(sim_number)
    """
    This function extracts the data from the saved jld2 files and places them in an excel.

    Arguments
    sim_number:                 The simulation number. It does not have to be the number, but it has to be how the map with data is saved.

    No returns, but The data that is extracted is:
    concentrations:                 A slice of all concentration through the middle of the simulation domain, at the time of saving
    bulk_concentrations:            The bulk concentrations at time of saving
    HRT:                            The HRT at the time of saving
    profiling:                      The time needed for each operation for each dT_bac
    nDiffIters:                     The amount of diffusion iteration required for each dT_bac step
    bulk_history:                   The bulk concentration after each dT_bac step
    Time:                           The time at each dT_bac step, dT_bac can be deducted from this
    """

    conc_saved, reactor_saved = load("results\\$(sim_number)\\results1D.jld2", "conc_saved", "reactor_saved")
    profiling, nDiffIters, bulk_history, Time = load("results\\$(sim_number)\\profilingResults.jld2", "profiling", "nDiffIters", "bulk_history", "Time")

    XLSX.openxlsx("results\\$(sim_number)\\julia_data_$(sim_number).xlsx", mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "Concentrations")
        for i in axes(conc_saved, 3)
            location = 1 + (i-1) * (size(conc_saved, 1) + 1)
            sheet["A$(location)"] = conc_saved[:,:,i]
        end

        sheet2 = XLSX.addsheet!(xf, "Bulk_concentrations")
        sheet2["A1"] = reactor_saved.bulk_concs

        sheet3 = XLSX.addsheet!(xf, "HRT")
        sheet3["A1"] = reactor_saved.HRT

        sheet4 = XLSX.addsheet!(xf, "Profiling")
        sheet4["B5"] = profiling

        sheet5 = XLSX.addsheet!(xf, "nDiffIters")
        sheet5["A1"] = nDiffIters

        sheet6 = XLSX.addsheet!(xf, "bulk_history")
        sheet6["A1"] = bulk_history

        sheet7 = XLSX.addsheet!(xf, "Time")
        sheet7["A1"] = Time.history
    end
    println("DONE!")
end
