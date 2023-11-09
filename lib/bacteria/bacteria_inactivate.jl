function bacteria_inactivate!(bac, constants)
    """
    This function inactivates any bacteria that are below the minimum mass threshold
    bac is a struct that contains all information regarding the bacteria
    constant is a struct containing all simulation constants    
    """
    mask_tooSmall = bac.molarMass * constants.bac_MW .< constants.min_bac_mass_grams # marks too small bacteria
    mask_positiveGrowthRate = bac.mu .> 0
    mask_possible_reactivation = .!bac.active .& mask_positiveGrowthRate # These bacteria were inactive but have a positive growth rate

    bac.active[mask_tooSmall .& bac.active] = 0 # active but too small --> inactivates
    # When possible 10% change of reactivation
    random_reactivation = rand(sum(mask_possible_reactivation)) .< 0.1 
    bac.active[mask_possible_reactivation] = random_reactivation

    # Reactivation is implicitely done already by inactive bacteria being able to grow.
    # If at any point those bacteria become larger than the threshold, they
    # will be set to active again.

    println("$(sum(random_reactivation)) individuals reactivated")

    return bac
end
