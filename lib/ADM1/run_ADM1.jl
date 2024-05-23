using DifferentialEquations
using ODE

function run_ADM1()  

    function adjust_parameter_for_temperature(base_parameter, enthalpy, R_joule, base_T, operational_T)
        adjusted_parameter = base_parameter * exp(enthalpy/R_joule * (1/base_T - 1/operational_T))
        return adjusted_parameter
    end

    # Operational parameters TODO
    Vliq = 1            # Liquid Volume         [m3]
    Vg = 1              # Gas volume            [m3]
    T = 308             # Temperature           [K]
    Pgas = 1            # Gas pressure          [bar]
    R = 8.314e-5        # Gas constant          [m3 bar/mol/K]
    R_joule = 8.314     # Gas constant          [J /mole/K]
    Q = 1               # Liquid inflowing      [m3/d]
    kla = 1             # Specific transfer rate d-1
    pH = 7.0
    S_H_ion = 10^-(pH) * 1000  # Proton concentration  [mol/m3]
    sim_time = 10       # Simulation duration  [d]

    # Rate parameters (specific for 35 Celsius)
    # Specific Monod Maximum uptake rates
    # UNITS NORMALIZED FOR THE RIGHT BIOMASS??
    # CHECK WHETHER THE 10x IS REQUIRED!!!!!!!!!!!
    km_sug = 250.0      # molS/molX/d
    km_but = 200.0      # molS/molX/d
    km_pro = 185.714    # molS/molX/d
    km_ac  = 200.0      # molS/molX/d
    km_h2  = 3200.0     # molS/molX/d

    #Monod half saturation constants (specific for 35 Celsius)
    K_S_sug = 2.604     # mol/m3
    K_S_but = 1.875     # mol/m3
    K_S_pro = 2.679     # mol/m3
    K_S_ac  = 2.34      # mol/m3
    K_S_h2  = 1.56e-3   # mol/m3

    # Decay rates (specific for 35 Celsius)
    kdec_sug = 0.2      # d-1
    kdec_but = 0.2      # d-1
    kdec_pro = 0.2      # d-1
    kdec_ac  = 0.2      # d-1
    kdec_h2  = 0.2      # d-1

    # Equilibrium constants TODO CHECK IF RIGHT TEMPERATURE!!!!!!!
    K_a_co2_base = 10^-6.35 * 1000  # [mol/m3]
    K_a_IN_base  = 10^-9.25 * 1000  # [mol/m3]
    K_a_ac  = 10^-4.76 * 1000       # [mol/m3]
    K_a_pr  = 10^-4.88 * 1000       # [mol/m3]
    K_a_bu  = 10^-4.82 * 1000       # [mol/m3]
    K_a_co2 = adjust_parameter_for_temperature(K_a_co2_base, 7646.0, R_joule, 298, T)
    K_a_IN =  adjust_parameter_for_temperature(K_a_IN_base, 51965.0, R_joule, 298, T)

    # UNITS NORMALIZED FOR THE RIGHT BIOMASS??   
    # Biomass yields (specific for 35 Celsius)
    Ysug = 0.12             # molX/molS
    Ybut = 0.06             # molX/molS
    Ypro = 0.028            # molX/molS
    Yac  = 0.02             # molX/molS
    Yh2  = 0.006            # molX/molS

    # Henry coefficients
    KH_h2_base  = 0.78       # mol/m3/bar_gas
    KH_ch4_base = 1.4        # mol/m3/bar_gas
    KH_co2_base = 35.0       # mol/m3/bar_gas
    KH_h2  = adjust_parameter_for_temperature(KH_h2_base,  -4180.0,  R_joule, 298, T)
    KH_ch4 = adjust_parameter_for_temperature(KH_ch4_base, -14240.0, R_joule, 298, T)
    KH_co2 = adjust_parameter_for_temperature(KH_co2_base, -19410.0, R_joule, 298, T)

    # Inhibition coefficients (specific for 35 Celsius)
    K_Ih2_but = 6.25e-4 # mol/m3
    K_Ih2_pro = 2.19e-4 # mol/m3
    K_I_nh3   = 1.80e-6 # mol/m3
    K_S_IN    = 1e-7    # mol/m3

    # Initial conditions (at t = 0) TODO
    S_sug_0 = 0         # mol/m3
    S_but_0 = 0         # mol/m3
    S_pro_0 = 0         # mol/m3
    S_ac_0  = 0         # mol/m3
    S_h2_0  = 0         # mol/m3
    S_ch4_0 = 0         # mol/m3
    S_IC_0  = 0         # mol/m3
    S_IN_0  = 0         # mol/m3

    X_sug_0 = 0         # mol/m3
    X_but_0 = 0         # mol/m3
    X_pro_0 = 0         # mol/m3
    X_ac_0  = 0         # mol/m3
    X_h2_0  = 0         # mol/m3
    X_I_0   = 0         # mol/m3

    S_h2_g_0  = 0       # mol/m3
    S_ch4_g_0 = 0       # mol/m3
    S_co2_g_0 = 0       # mol/m3

    c0 = [S_sug_0, S_but_0, S_pro_0, S_ac_0, S_h2_0, S_ch4_0, S_IC_0, S_IN_0,
            X_sug_0, X_but_0, X_pro_0, X_ac_0, X_h2_0, X_I_0,
            S_h2_g_0, S_ch4_g_0, S_co2_g_0]

            # CHECK WHETHER THESE NUMBERS ARE STILL CORRECT!!!!!!!!!!

                    #   Gluc,    But,               Pro,                     Ac,                H2,              CH4,          IC,                      IN,      Xsug,   Xbut,  Xpro,  Xac,   Xh2,   XI
    reaction_matrix = [-1        0.156-0.13*Ysug    0.462857-0.38571*Ysug    1.23-1.025*Ysug    2.28-1.9*Ysug    0             1.527429-1.27286*Ysug   -Ysug     Ysug    0      0      0      0      0;   # Glucose uptake
                        0       -1                  0                        2-2.5*Ybut         2                0             0                       -Ybut     0       Ybut   0      0      0      0;   # Butyrate uptake
                        0        0                 -1                        1-1.67*Ypro        3-3.33*Ypro      0             1-1.67*Ypro             -Ypro     0       0      Ypro   0      0      0;   # Propionate uptake
                        0        0                  0                       -1                  0                1-2.5*Yac     1-2.5*Yac               -Yac      0       0      0      Yac    0      0;   # Acetate uptake
                        0        0                  0                        0                 -1                0.25-2.5*Yh2 -0.25-2.5*Yh2            -Yh2      0       0      0      0      Yh2    0;   # H2 uptake
                        0        0                  0                        0                  0                0             0                        0       -1       0      0      0      0      1;   # X_sug decay
                        0        0                  0                        0                  0                0             0                        0        0      -1      0      0      0      1;   # X_but decay
                        0        0                  0                        0                  0                0             0                        0        0       0     -1      0      0      1;   # X_pro decay
                        0        0                  0                        0                  0                0             0                        0        0       0      0     -1      0      1;   # X_ac  decay
                        0        0                  0                        0                  0                0             0                        0        0       0      0      0     -1      1]   # X_h2  decay

    function determine_inhibition(pH, UL, LL)
        if pH > UL
            return 1.0
        else
            I = exp(-3 * ((pH - UL) / (UL-LL))^2)
            return I
        end
    end

    function ADM1_rates(concs)

        S_sug, S_but_tot, S_pro_tot, S_ac_tot, S_h2, S_ch4, S_IC, S_IN,
            X_sug, X_but, X_pro, X_ac, X_h2, X_I,
            S_h2_g, S_ch4_g, S_co2_g = concs


        I_pH_aa = determine_inhibition(pH, 5.5, 4.0)    # [-]
        I_pH_ac = determine_inhibition(pH, 7.0, 6.0)    # [-]
        I_pH_h2 = determine_inhibition(pH, 6.0, 5.0)    # [-]

        S_nh3 = S_IN * K_a_IN / (S_H_ion + K_a_IN)      # [mol/m3]

        # Inhibitors
        I_IN_lim = 1 / (1 + K_S_IN/S_IN)            # [-]
        I_h2_but = 1 / (1 + S_h2/K_Ih2_but)         # [-]
        I_h2_pro = 1 / (1 + S_h2/K_Ih2_pro)         # [-]
        I_nh3    = 1 / (1 + S_nh3/K_I_nh3)          # [-]

        I1 = I_pH_aa * I_IN_lim             # [-]
        I2 = I1 * I_h2_but                  # [-]
        I3 = I1 * I_h2_pro                  # [-]
        I4 = I_pH_ac * I_IN_lim * I_nh3     # [-]
        I5 = I_pH_h2 * I_IN_lim             # [-]

        S_but = S_but_tot * S_H_ion / (K_a_bu + S_H_ion)    # [mol/m3]
        S_pro = S_pro_tot * S_H_ion / (K_a_pr + S_H_ion)    # [mol/m3]
        S_ac  = S_ac_tot  * S_H_ion / (K_a_ac + S_H_ion)    # [mol/m3]

        rho1  = km_sug * S_sug / (K_S_sug + S_sug) * X_sug * I1         # sugar uptake           [mol/m3/d]
        rho2  = km_but * S_but / (K_S_but + S_but) * X_but * I2         # butyrate uptake        [mol/m3/d]
        rho3  = km_pro * S_pro / (K_S_pro + S_pro) * X_pro * I3         # propionate uptake      [mol/m3/d]
        rho4  = km_ac  * S_ac  / (K_S_ac  + S_ac ) * X_ac  * I4         # acetate uptake         [mol/m3/d]
        rho5  = km_h2  * S_h2  / (K_S_h2  + S_h2 ) * X_h2  * I5         # hydrogen uptake        [mol/m3/d]
        rho6  = kdec_sug * X_sug            # X_sug Decay   [mol/m3/d]
        rho7  = kdec_but * X_but            # X_but Decay   [mol/m3/d]
        rho8  = kdec_pro * X_pro            # X_pro Decay   [mol/m3/d]
        rho9  = kdec_ac  * X_ac             # X_ac  Decay   [mol/m3/d]
        rho10 = kdec_h2  * X_h2             # X_h2  Decay   [mol/m3/d]

        rho_array = [rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8, rho9, rho10]
        liquid_rates = dropdims(sum(reaction_matrix .* rho_array, dims=1), dims=1)  # [moli/m3/d]
        return liquid_rates
    end
    
    function ADM1_balances(concs, p, t)

        S_sug, S_but, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN,
            X_sug, X_but, X_pro, X_ac, X_h2, X_I,
            S_h2_g, S_ch4_g, S_co2_g = concs

        S_co2 = S_IC * S_H_ion / (S_H_ion + K_a_co2)    # [mol/m3]
        p_h2o = 0.0313 * exp(43980/8.314 * (1/298 - 1/T)) # [bar]

        # Liquid rates
        liquid_rates = ADM1_rates(concs)   # [moli/m3/d]

        # Solubilities
        S_h2_s  = KH_h2  * S_h2_g  * R * T      # [mol/m3]
        S_ch4_s = KH_ch4 * S_ch4_g * R * T      # [mol/m3]
        S_co2_s = KH_co2 * S_co2_g * R * T      # [mol/m3]

        # Gas transfer rates
        rt_h2  = kla * (S_h2  - S_h2_s )        # [mol/m3/d]
        rt_ch4 = kla * (S_ch4 - S_ch4_s)        # [mol/m3/d]
        rt_co2 = kla * (S_co2 - S_co2_s)        # [mol/m3/d]
        Qg = R * T / (Pgas - p_h2o) * Vliq * (rt_h2 + rt_ch4 + rt_co2)  # [m3/d]          # Gas flow m3/d

        # Inflow concentrations TODO
        # S_sug_in, S_but_in, S_pro_in, S_ac_in, S_H2_in, S_ch4_in, S_co2_in, S_nh3_in, X_sug_in, X_but_in, X_pro_in, X_ac_in, X_H2_in, X_I_in
        conc_in = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]    # [mol/m3]

        # SOLVE pH???? Does not seem necessary as we enfore a pH of 7 (and everything is bulk)
        # Use similar solve system as the IbM (control_pH in calculate_bulk_concentrations)

        # All liquid compound balances
        f = zeros(size(concs))
        f[1:length(liquid_rates)] = Q ./ Vliq .* (conc_in .- concs[1:length(liquid_rates)]) .+ liquid_rates # [mol/m3/d]

        # Adjustment for gas transfer
        f[5] = f[5] - rt_h2     # [mol/m3/d]
        f[6] = f[6] - rt_ch4    # [mol/m3/d]
        f[7] = f[7] - rt_co2    # [mol/m3/d]

        # Gas compound balances
        f[15] = - S_h2_g  * Qg / Vg + rt_h2  * Vliq/Vg  # [mol/m3/d]
        f[16] = - S_ch4_g * Qg / Vg + rt_ch4 * Vliq/Vg  # [mol/m3/d]
        f[17] = - S_co2_g * Qg / Vg + rt_co2 * Vliq/Vg  # [mol/m3/d]

        return f
    end

    prob = ODEProblem(ADM1_balances, c0, (0.0, sim_time))
    sol = solve(prob, Rosenbrock23(autodiff=false), isoutofdomain=(y,p,t)->any(x->x.<0,y), reltol=1e-8, abstol=1e-20)
    return sol.u
end
