function bacteria_inactivate!(bac_vecfloat, bac_vecbool, constants)
    """
    This function inactivates any bacteria that are below the minimum mass threshold
    When bacteria are not active but have positive growth, they have a 10 % change to become active again
    
    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    constants:          A "General" struct containing all the simulation constants  

    Returns
    bac:                A bac struct with updates activation status
    """
    mask_tooSmall = bac_vecfloat.molarMass * constants.bac_MW .< constants.min_bac_mass_grams    # marks too small bacteria
    mask_positiveGrowthRate = bac_vecfloat.mu .> 0
    mask_possible_reactivation = .!bac_vecbool.active .& mask_positiveGrowthRate                # These bacteria were inactive but have a positive growth rate

    bac_vecbool.active[mask_tooSmall .& bac_vecbool.active] .= 0                                        # active but too small --> inactivates

    # When possible 10% change of reactivation
    activation_chance = rand(sum(mask_possible_reactivation))
    random_reactivation = activation_chance .< 0.1 
    bac_vecbool.active[mask_possible_reactivation] = random_reactivation

    # Reactivation is implicitely done already by inactive bacteria being able to grow.
    # If at any point those bacteria become larger than the threshold, they
    # will be set to active again.

    println("$(sum(random_reactivation)) individuals reactivated")

    return bac_vecbool
end
