function calculate_inhibition(UL, LL, pH)
    """
    This function calculates the inhibition factor based on the pH and an empirical equation.

    Arguments
    UL:                     The Upper Limit, at this pH there is no inhibition
    LL:                     The Lower Limit, at this pH there is total inhibition
    pH:                     The pH of the grid cell

    Returns
    inhibition:             The inhibition factor
    """

    if pH > UL
        inhibition = 1
    else
        inhibition = exp(-3 * ((pH - UL)/(UL - LL)) ^2)
    end

    return inhibition
end

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

    pH = -log10(Sh)

    if species == 1 # BO
        UL = 5.5
        LL = 4.0
        maint = 0.2/24 # h-1
        mu_max = 1.2 * calculate_inhibition(UL, LL, pH) / 24

    elseif species == 2 # AM
        UL = 7.0
        LL = 6.0
        maint = 0.2/24 # h-1
        mu_max = 0.4 * calculate_inhibition(UL, LL, pH) / 24

        # AM Alternative, upper and lower inhibition
        # UL = 9.17
        # LL = 5.95
        # mu_max = 0.4 * (1+2*10^(0.5*(LL-UL)))/(1+10^(pH-UL)+10^(LL-pH)) /24
    
    elseif species == 3 # HM
        UL = 6.0
        LL = 5.0
        maint = 0.2/24 # h-1
        mu_max = 2.1 * calculate_inhibition(UL, LL, pH) / 24
    
    # Nitrospira values
    # if species == 1 # AOB
    #     mu_max = ((1.28*10^(12) * exp(-8183/T)) / (1 + ((2.05*10^(-9))/Sh) + (Sh/(1.66*10^(-7)))))/24
    #     maint = (1.651*10^(11) * exp(-8183/T))/24
        
    # # elseif species == 2 # NOB Nitrobacter
    # #     mu_max = ((6.69*10^(7) * exp(-5295/T)) / (1 + ((2.05*10^(-9))/Sh) + (Sh/(1.66*10^(-7)))))/24
    # #     maint = (8.626*10^(6)*exp(-5295/T))/24

    # elseif species == 2 # NOB Nitrospira
    #     mu_max = 0.63 * ((6.69*10^(7) * exp(-5295/T)) / (1 + ((2.05*10^(-9))/Sh) + (Sh/(1.66*10^(-7)))))/24
    #     maint = 0.63 * (8.626*10^(6)*exp(-5295/T))/24

    # elseif species == 3 # AMX Brocadia spp [Puyol 2014]
    #     mu_max = 1.89*10^(8) * exp(-7330/T)
    #     maint = 0.05 * mu_max

    else
        error("Bacterial species not implemented: $(species)")
    end

    return mu_max, maint    
end
