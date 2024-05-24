def jorek_to_si(central_density, central_mass, name, jor_value, shownames):

    import math
    
    # Constants
    pi             = 3.14159265359;
    ech            = 1.60217662e-19;
    mu0            = 4.0*pi*1.0e-7;
    mproton        = 1.6726219e-27;
    n0             = central_density * 1e20;
    rho0           = central_mass * n0 * mproton;
    sqrt_rhomu     = math.sqrt(rho0*mu0);
    sqrt_mu_ov_rho = math.sqrt(mu0/rho0);
    
    # Cases for each parameter
    
    if (name == 'eta'):
        si_val = jor_value * sqrt_mu_ov_rho

    if (name == 'T_eV'):
        si_val = jor_value / (ech * mu0 * n0)

    if (name == 'pressure'):
        si_val = jor_value / mu0

    if (name == 'velocity'):
        si_val = jor_value / sqrt_rhomu

    if (name == 'time'):
        si_val = jor_value * sqrt_rhomu

    if (name == 'gr_rate'):
        si_val = jor_value / sqrt_rhomu
        
    if (name == 'Efield'):
        si_val = jor_value * sqrt_rhomu

    if (name == 'D_perp'):
        si_val = jor_value / sqrt_rhomu

    if (name == 'visco_dyn'):
        si_val = jor_value / sqrt_mu_ov_rho
        
    if (name == 'visco_kin'):
        si_val = jor_value / (sqrt_mu_ov_rho * rho0)
        
    if (name == 'kappa_dyn'):
        si_val = jor_value / sqrt_mu_ov_rho
        print("You still need to multiply by (gamma-1)")

    if (name == 'kappa_m2/s'):
        si_val = jor_value / (sqrt_mu_ov_rho * rho0)
        print("You still need to multiply by (gamma-1)")

    if (name == 'kappa_ms-1'):
        si_val = jor_value / (sqrt_mu_ov_rho * rho0 / n0)
        print("You still need to multiply by (gamma-1)")
        
        
    if (shownames):
        print("eta")
        print("T_eV")
        print("pressure")
        print("velocity")
        print("time")
        print("gr_rate")
        print("Efield")
        print("D_perp")
        print("kappa_dyn")
        print("kappa_m2/s")
        print("kappa_ms-1")
        print("visco_dyn")
        print("visco_kin")
    
        
    return si_val



def to_jorek_units(central_density, central_mass, name, si_value, shownames):

    import math
    
    # Constants
    pi             = 3.14159265359;
    ech            = 1.60217662e-19;
    mu0            = 4.0*pi*1.0e-7;
    mproton        = 1.6726219e-27;
    n0             = central_density * 1e20;
    rho0           = central_mass * n0 * mproton;
    sqrt_rhomu     = math.sqrt(rho0*mu0);
    sqrt_mu_ov_rho = math.sqrt(mu0/rho0);
    
    # Cases for each parameter
    
    if (name == 'eta'):
        jor_val = si_value / sqrt_mu_ov_rho

    if (name == 'T_eV'):
        jor_val = si_value * ech * mu0 * n0

    if (name == 'pressure'):
        jor_val = si_value * mu0

    if (name == 'velocity'):
        jor_val = si_value * sqrt_rhomu

    if (name == 'time'):
        jor_val = si_value / sqrt_rhomu

    if (name == 'gr_rate'):
        jor_val = si_value * sqrt_rhomu
        
    if (name == 'Efield'):
        jor_val = si_value / sqrt_rhomu

    if (name == 'D_perp'):
        jor_val = si_value * sqrt_rhomu

    if (name == 'visco_dyn'):
        jor_val = si_value * sqrt_mu_ov_rho
        
    if (name == 'visco_kin'):
        jor_val = si_value * sqrt_mu_ov_rho * rho0
        
    if (name == 'kappa_dyn'):
        jor_val = si_value * sqrt_mu_ov_rho
        print("You still need to multiply by (gamma-1)")

    if (name == 'kappa_m2/s'):
        jor_val = si_value * sqrt_mu_ov_rho * rho0
        print("You still need to multiply by (gamma-1)")

    if (name == 'kappa_ms-1'):
        jor_val = si_value * sqrt_mu_ov_rho * rho0 / n0
        print("You still need to multiply by (gamma-1)")
        
        
    if (shownames):
        print("eta")
        print("T_eV")
        print("pressure")
        print("velocity")
        print("time")
        print("gr_rate")
        print("Efield")
        print("D_perp")
        print("kappa_dyn")
        print("kappa_m2/s")
        print("kappa_ms-1")
        print("visco_dyn")
        print("visco_kin")
        
    return jor_val