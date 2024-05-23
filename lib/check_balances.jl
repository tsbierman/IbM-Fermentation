function check_balances(bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint,constants_vecstring, constants_vecbool, constants_matfloat, settings_string, reaction_matrix, bulk_concentrations, invHRT, bulk_change, tolerance)
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

    estimated_compounds = zeros(length(constants_vecstring.compoundNames[constants_vecint.Gas_k .!= 1]))
    Sh = 10^-constants_float.pHsetpoint

    for specie in 1:length(constants_vecstring.speciesNames)
        specie_index = (bac_vecint.species .== specie) .* bac_vecbool.active                    # Find index of active bacteria of specific specie
        specie_mu_withMain = bac_vecfloat.mu[specie_index]                                       # Obtain their growth rates
        
        # Decay change
        decay_index = specie_mu_withMain .< 0                                           # Find which ones are decaying but still active
        mass_decaying_bac = bac_vecfloat.molarMass[specie_index][decay_index]                    # Select mass [mol] of decaying bacteria of the specie
        decayed_biomass = sum(mass_decaying_bac .* specie_mu_withMain[decay_index])     # Calculate how much biomass in total [mol_X] is lost due to decay
        compoundChange_decay = constants_matfloat.MatrixDecay[:, specie] .* decayed_biomass       # Compound change [mol_i/h], note: NEGATIVE!

        # Growth change
        if isnan(constants_vecfloat.maintenance[1])
            _, specie_maint = determine_max_growth_rate_and_maint(specie, constants_float.T, Sh)
        else
            specie_maint = constants_vecfloat.maintenance[specie]                                    # Get maintenance
        end
        
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

    CO2_index = findall(constants_vecstring.compoundNames.== "CO2")[1]
    Keq = constants_matfloat.Keq
    Denm = (1 .+ Keq[CO2_index, 1]) .* Sh^3 .+ Keq[CO2_index, 2] .* Sh^2 .+ Keq[CO2_index, 2] .* Keq[CO2_index, 3] .* Sh .+ Keq[CO2_index, 2] .* Keq[CO2_index, 3] .* Keq[CO2_index, 4]
    CO2_conc = (bulk_concentrations[CO2_index] .* Sh^3) ./ Denm  

    copybulk_conc = copy(outflow)
    copybulk_conc[CO2_index] = CO2_conc

    # Calculations required for the gas-liquid transfer
    p_h2o = 0.0313 * exp(43980/(constants_float.R*1000) * (1/298 - 1/constants_float.T))
    solubilities = constants_vecfloat.Kh[constants_vecint.Gas_k .== 1] .* outflow[constants_vecint.Gas_k .== 1] .* constants_float.R*10 .* constants_float.T       # Same units as prev_conc (as long as R and Kh are in the same units [mol/L]
    gas_transfer_rates = constants_float.kla .* (copybulk_conc[constants_vecint.Gas_k .== -1] .- solubilities) # [mol/L/h]
    Qgas = constants_float.R*10 * constants_float.T / (constants_float.Pgas - p_h2o) * constants_float.Vr * sum(gas_transfer_rates) 

    # Conversion factor volume slice to volume sphere
    f = calculate_slice_sphere_conversion(bac_vecfloat, bac_vecbool, constants_float, settings_string)
    actual_compoundChange = actual_compoundChange .* f ./ constants_float.Vr    # Convert from [mol_i/h] to [mol_i/L/h] (reactor level)
    actual_compoundChange = [actual_compoundChange; zeros(sum(constants_vecint.Gas_k .== 1))] # Patch for gaseous compounds

    balances = zeros(size(outflow))

    change_liquid = .!constants_vecbool.Dir_k .& (constants_vecint.Gas_k .!= 1)

    balances[change_liquid] = (influent[change_liquid] .- outflow[change_liquid]) .* invHRT .+ actual_compoundChange[change_liquid] .- bulk_change[change_liquid]
    balances[constants_vecint.Gas_k .== -1] = balances[constants_vecint.Gas_k .== -1] .- gas_transfer_rates
    balances[constants_vecint.Gas_k .== 1] = -Qgas ./ constants_float.Vgas .* outflow[constants_vecint.Gas_k .== 1] .+ gas_transfer_rates .* constants_float.Vr ./ constants_float.Vgas .- bulk_change[constants_vecint.Gas_k .== 1]

    balances_closes = maximum(abs.(balances[.!constants_vecbool.Dir_k])) <= tolerance
    dirichlet_closes = maximum(abs.(balances[constants_vecbool.Dir_k] .- actual_compoundChange[constants_vecbool.Dir_k])) <= tolerance

    return biomass_closes, balances_closes, dirichlet_closes
end
