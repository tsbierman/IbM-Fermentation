using DifferentialEquations
using ODE
using Plots

function run_ADM1()  

    function adjust_parameter_for_temperature(base_parameter, enthalpy, R_joule, base_T, operational_T)
        adjusted_parameter = base_parameter * exp(enthalpy/R_joule * (1/base_T - 1/operational_T))
        return adjusted_parameter
    end

    # Operational parameters
    Vliq = 45           # Liquid Volume          [L]
    Vg = 15             # Gas volume             [L]
    T = 308.15          # Temperature            [K]
    Pgas = 1.013        # Gas pressure           [bar]
    R = 8.314e-2        # Gas constant           [L bar/mol/K]
    R_joule = 8.314     # Gas constant           [J/mole/K]
    Q = 0.1875          # Liquid inflowing       [L/h] Can adjust this number to play with HRT
    kla = 200/24        # Specific transfer rate [h-1]
    pH = 7.0
    S_H_ion = 10^-(pH)  # Proton concentration   [mol/L]
    sim_time = 50*24    # Simulation duration    [h]

    # Rate parameters (specific for 35 Celsius)
    # Specific Monod Maximum uptake rates
    km_but = 4/24      # molS/molX/h
    km_ac  = 4/24      # molS/molX/h
    km_h2  = 70/24     # molS/molX/h

    #Monod half saturation constants (specific for 35 Celsius)
    K_S_but = 1.8750e-3     # mol/L
    K_S_ac  = 2.3438e-3     # mol/L
    K_S_h2  = 1.5625e-6     # mol/L

    # Decay rates: 10% of the maximum growth rate
    kdec_but = 1.2/10/24      # h-1
    kdec_ac  = 0.4/10/24      # h-1
    kdec_h2  = 2.1/10/24      # h-1

    # Equilibrium constants
    K_a_co2_base = 10^-6.35  # [mol/L]
    K_a_IN_base  = 10^-9.25  # [mol/L]
    K_a_ac  =      10^-4.76  # [mol/L]
    K_a_bu  =      10^-4.82  # [mol/L]
    K_a_co2 = adjust_parameter_for_temperature(K_a_co2_base, 7646.0, R_joule, 298, T)
    K_a_IN =  adjust_parameter_for_temperature(K_a_IN_base, 51965.0, R_joule, 298, T)

    # Biomass yields (specific for 35 Celsius)
    Ybut = 0.3             # molX/molS
    Yac  = 0.1             # molX/molS
    Yh2  = 0.03            # molX/molS

    # Henry coefficients
    KH_h2_base  = 0.00078       # mol/L/bar_gas
    KH_ch4_base = 0.0014        # mol/L/bar_gas
    KH_co2_base = 0.035         # mol/L/bar_gas
    KH_h2  = adjust_parameter_for_temperature(KH_h2_base,  -4180.0,  R_joule, 298, T)
    KH_ch4 = adjust_parameter_for_temperature(KH_ch4_base, -14240.0, R_joule, 298, T)
    KH_co2 = adjust_parameter_for_temperature(KH_co2_base, -19410.0, R_joule, 298, T)

    # Inhibition coefficients (specific for 35 Celsius)
    K_Ih2_but = 6.25e-7 # mol/L
    K_I_nh3   = 1.80e-3 # mol/L

    # Initial conditions (at t = 0)
    S_but_0 = 0.05      # mol/L
    S_ac_0  = 0         # mol/L
    S_h2_0  = 0         # mol/L
    S_ch4_0 = 0         # mol/L
    S_IC_0  = 0.001     # mol/L
    S_IN_0  = 0.03611   # mol/L

    X_but_0 = 0.368         # mol/L  25 g/L total biomass
    X_ac_0  = 0.368         # mol/L 
    X_h2_0  = 0.368         # mol/L 
    X_I_0   = 0             # mol/L

    S_h2_g_0  = 0       # mol/L
    S_ch4_g_0 = 0       # mol/L
    S_co2_g_0 = 0       # mol/L

    c0 = [S_but_0, S_ac_0, S_h2_0, S_ch4_0, S_IC_0, S_IN_0,
            X_but_0, X_ac_0, X_h2_0, X_I_0,
            S_h2_g_0, S_ch4_g_0, S_co2_g_0]

                    #   But,     Ac,          H2,       CH4,             IC,            IN,          Xbut,  Xac,   Xh2,   XI
    reaction_matrix = [-1        2-0.5*Ybut   2         0                0             -0.2*Ybut     Ybut   0      0      0;   # Butyrate uptake
                        0       -1            0         1-0.5*Yac        1-0.5*Yac     -0.2*Yac      0      Yac    0      0;   # Acetate uptake
                        0        0           -1         0.25-0.5*Yh2    -0.25-0.5*Yh2  -0.2*Yh2      0      0      Yh2    0;   # H2 uptake
                        0        0            0         0                0              0           -1      0      0      1;   # X_but decay
                        0        0            0         0                0              0            0     -1      0      1;   # X_ac  decay
                        0        0            0         0                0              0            0      0     -1      1]   # X_h2  decay

    function determine_inhibition(pH, UL, LL)
        if pH > UL
            return 1.0
        else
            I = exp(-3 * ((pH - UL) / (UL-LL))^2)
            return I
        end
    end

    function ADM1_rates(concs)

        S_but_tot, S_ac_tot, S_h2, S_ch4, S_IC, S_IN,
            X_but, X_ac, X_h2, X_I,
            S_h2_g, S_ch4_g, S_co2_g = concs

        I_pH_aa = determine_inhibition(pH, 5.5, 4.0)    # [-]
        I_pH_ac = determine_inhibition(pH, 7.0, 6.0)    # [-]
        I_pH_h2 = determine_inhibition(pH, 6.0, 5.0)    # [-]

        S_nh3 = S_IN * K_a_IN / (S_H_ion + K_a_IN)      # [mol/L]

        # Inhibitors
        I_h2_but = 1 / (1 + S_h2/K_Ih2_but)         # [-]
        I_nh3    = 1 / (1 + S_nh3/K_I_nh3)          # [-]

        I1 = I_pH_aa * I_h2_but  # [-]
        I2 = I_pH_ac * I_nh3     # [-]
        I3 = I_pH_h2             # [-]

        S_but = S_but_tot * K_a_bu / (K_a_bu + S_H_ion)    # [mol/L] Uses But-
        S_ac  = S_ac_tot  * K_a_ac / (K_a_ac + S_H_ion)    # [mol/L] Uses Ac-

        rho1 = km_but * S_but / (K_S_but + S_but) * X_but * I1         # butyrate uptake        [mol/L/h]
        rho2 = km_ac  * S_ac  / (K_S_ac  + S_ac ) * X_ac  * I2         # acetate uptake         [mol/L/h]
        rho3 = km_h2  * S_h2  / (K_S_h2  + S_h2 ) * X_h2  * I3         # hydrogen uptake        [mol/L/h]
        rho4 = kdec_but * X_but                                        # X_but Decay   [mol/L/h]
        rho5 = kdec_ac  * X_ac                                         # X_ac  Decay   [mol/L/h]
        rho6 = kdec_h2  * X_h2                                         # X_h2  Decay   [mol/L/h]

        rho_array = [rho1, rho2, rho3, rho4, rho5, rho6]
        liquid_rates = dropdims(sum(reaction_matrix .* rho_array, dims=1), dims=1)  # [moli/L/h]
        return liquid_rates
    end
    
    function ADM1_balances(concs, p, t)

        S_but, S_ac, S_h2, S_ch4, S_IC, S_IN,
            X_but, X_ac, X_h2, X_I,
            S_h2_g, S_ch4_g, S_co2_g = concs

        S_co2 = S_IC * S_H_ion / (S_H_ion + K_a_co2)      # [mol/L]
        p_h2o = 0.0313 * exp(43980/8.314 * (1/298 - 1/T)) # [bar]

        # Liquid rates
        liquid_rates = ADM1_rates(concs)   # [moli/L/h]

        # Solubilities
        S_h2_s  = KH_h2  * S_h2_g  * R * T      # [mol/L]
        S_ch4_s = KH_ch4 * S_ch4_g * R * T      # [mol/L]
        S_co2_s = KH_co2 * S_co2_g * R * T      # [mol/L]

        # Gas transfer rates
        rt_h2  = kla * (S_h2  - S_h2_s )        # [mol/L/h]
        rt_ch4 = kla * (S_ch4 - S_ch4_s)        # [mol/L/h]
        rt_co2 = kla * (S_co2 - S_co2_s)        # [mol/L/h]
        Qg = R * T / (Pgas - p_h2o) * Vliq * (rt_h2 + rt_ch4 + rt_co2)  # [L/h]

        # Inflow concentrations
        # S_but_in, S_ac_in, S_H2_in, S_ch4_in, S_co2_in, S_nh3_in, X_but_in, X_ac_in, X_H2_in, X_I_in
        conc_in = [0.05, 0, 0, 0, 0.001, 0.03611, 0, 0, 0, 0]    # [mol/L]

        # All liquid compound balances
        f = zeros(size(concs))
        f[1:length(liquid_rates)] = Q ./ Vliq .* (conc_in .- concs[1:length(liquid_rates)]) .+ liquid_rates # [mol/L/h]

        # Adjustment for gas transfer
        f[3] = f[3] - rt_h2     # [mol/L/h]
        f[4] = f[4] - rt_ch4    # [mol/L/h]
        f[5] = f[5] - rt_co2    # [mol/L/h]

        # Gas compound balances
        f[11] = - S_h2_g  * Qg / Vg + rt_h2  * Vliq/Vg  # [mol/L/h]
        f[12] = - S_ch4_g * Qg / Vg + rt_ch4 * Vliq/Vg  # [mol/L/h]
        f[13] = - S_co2_g * Qg / Vg + rt_co2 * Vliq/Vg  # [mol/L/h]

        return f
    end

    prob = ODEProblem(ADM1_balances, c0, (0.0, sim_time))
    sol = solve(prob, Rodas4P(autodiff=false), isoutofdomain=(y,p,t)->any(x->x.<0,y), reltol=1e-8, abstol=1e-20)
    return sol.t, sol.u
end
