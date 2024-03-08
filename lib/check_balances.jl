function check_balances(bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecstring, constants_vecbool, constants_matfloat, settings_string, reaction_matrix, bulk_concentrations, invHRT, bulk_change, tolerance)
    """
    This function checks whether the amount of biomass produces/decayed alligns with the compounds consumed/produced.
    It also checks whether the balance of the whole reactor closes.

    Arguments
    bac:                    A "General" struct containing all parameters related to the bacteria
    constants:              A "General" struct containing all the simulation constants
    settings:               A "General" struct containing all the settings of the simulation
    reaction_matrix:        A (ny, nx, ncompounds) matrix containing all reaction rates per gridcell and compound [mol/L/h]
    bulk_concentrations:    A (ncompounds,) Vector that contains the bulk concentrations [mol/L]
    invHRT:                 The inverse of the HRT, equal to the dilution rate
    bulk_change:            The change of the bulk concentration compared to before [mol/L/h]
    tolerance:              The tolerance below which the residuals of the balances should be

    Returns
    biomass_closes:         A Boolean indicating whether the compounds produced due to growth/decay allign with what is expected
    balance_closes:         A Boolean indicating whether the overall reactor balances close
    dirichlet_closes:       A Boolean indicating whether the dirichlet boundary condition balances close
    """

    estimated_compounds = zeros(length(constants_vecstring.compoundNames))
    
    for specie in 1:length(constants_vecstring.speciesNames)
        specie_index = (bac_vecint.species .== specie) .* bac_vecbool.active                    # Find index of active bacteria of specific specie
        specie_mu_withMain = bac_vecfloat.mu[specie_index]                                       # Obtain their growth rates
        
        # Decay change
        decay_index = specie_mu_withMain .< 0                                           # Find which ones are decaying but still active
        mass_decaying_bac = bac_vecfloat.molarMass[specie_index][decay_index]                    # Select mass [mol] of decaying bacteria of the specie
        decayed_biomass = sum(mass_decaying_bac .* specie_mu_withMain[decay_index])     # Calculate how much biomass in total [mol_X] is lost due to decay
        compoundChange_decay = constants_matfloat.MatrixDecay[:, specie] .* decayed_biomass       # Compound change [mol_i/h], note: NEGATIVE!

        # Growth change
        specie_maint = constants_vecfloat.maintenance[specie]                                    # Get maintenance
        specie_mu_withoutMain = specie_mu_withMain .+ specie_maint                      # Convert mu
        increased_biomass = sum(bac_vecfloat.molarMass[specie_index] .* specie_mu_withoutMain)   # Increase of biomass due to growth
        compoundChange_growth = constants_matfloat.MatrixMet[:,specie] .* increased_biomass       # Compound change due to growth   [mol_i/h]

        # Update estimate
        estimated_compounds = estimated_compounds .+ compoundChange_growth .- compoundChange_decay # Update estimated compounds, - compound_decay due to it begin negative        
    end

    actual_compoundChange = dropdims(sum(reaction_matrix, dims=(1,2)), dims=(1,2)) .* constants_float.Vg       # Convert from [mol_i/L/h] (grid level) to [mol_i/h]
    difference = actual_compoundChange .- estimated_compounds
    biomass_closes = maximum(abs.(difference)) <= tolerance

    influent = constants_vecfloat.influent_concentrations
    outflow = bulk_concentrations

    # Conversion factor volume slice to volume sphere
    f = calculate_slice_sphere_conversion(bac_vecfloat, bac_vecbool, constants_float, settings_string)
    actual_compoundChange = actual_compoundChange .* f ./ constants_float.Vr    # Convert from [mol_i/h] to [mol_i/L/h] (reactor level)

    balances = (influent .- outflow) .* invHRT .+ actual_compoundChange .- bulk_change
    balances_closes = maximum(abs.(balances[.!constants_vecbool.Dir_k])) <= tolerance
    dirichlet_closes = maximum(abs.(balances[constants_vecbool.Dir_k] .- actual_compoundChange[constants_vecbool.Dir_k])) <= tolerance

    return biomass_closes, balances_closes, dirichlet_closes
end
