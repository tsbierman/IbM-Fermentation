function determine_max_growth_rate_and_maint(species, T, Sh)
    """
    This function determines the maximum growth rate and maintenance [h-1] for a 
    specific species under certain conditions
    Take note that this function is completely hardcoded
    This function is only called upon when there are no mu_max and maintenance supplied in the excel

    Arguments
    species:                An integer that indicates the bacteria species
    T:                      The temperature (Kelvin)
    Sh:                     The concentration of protons

    Returns
    mu_max:                 The calculated maximum growth rate (mu) for this specie
    maint:                  The calculated maintenance for this specie
    """

    if species == 1 # AOB
        mu_max = ((1.28*10^(12) * exp(-8183/T)) / (1 + ((2.05*10^(-9))/Sh) + (Sh/(1.66*10^(-7)))))/24
        maint = (1.651*10^(11) * exp(-8183/T))/24
        
    # elseif species == 2 # NOB Nitrobacter
    #     mu_max = ((6.69*10^(7) * exp(-5295/T)) / (1 + ((2.05*10^(-9))/Sh) + (Sh/(1.66*10^(-7)))))/24
    #     maint = (8.626*10^(6)*exp(-5295/T))/24

    elseif species == 2 # NOB Nitrospira
        mu_max = 0.63 * ((6.69*10^(7) * exp(-5295/T)) / (1 + ((2.05*10^(-9))/Sh) + (Sh/(1.66*10^(-7)))))/24
        maint = 0.63 * (8.626*10^(6)*exp(-5295/T))/24

    elseif species == 3 # AMX Brocadia spp [Puyol 2014]
        mu_max = 1.89*10^(8) * exp(-7330/T)
        maint = 0.05 * mu_max

    else
        error("Bacterial species not implemented: $(species)")
    end

    return mu_max, maint    
end
