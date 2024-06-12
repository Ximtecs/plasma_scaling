include("scaling_CGS.jl")

#----------- strcut for code units needed to reconstuct real units ----
struct code_units
    c::Float64
    B::Float64
    p::Float64
    rho::Float64
    debye::Float64
end
#------------------------------------------------

#-------- struct with the single real value needed to reconstruct all real units ----
struct known_real_units
    number_density ::Float64
end
#------------------------------------------------


#-------- strcut for fudging physical valeus -----------
struct ratios_and_scaling
    electron_proton_ratio :: Float64
    eps_0_scaling :: Float64
    q_charge_scaling :: Float64
end
#------------------------------------------------


function find_real_units_CGS(code_units, known_real_units, ratios_and_scaling)
    
    base = BaseUnits("CGS")
    #------------- Scale base units ----------------
        #------ calculate electron mass scaing based on ratio of electron to proton mass -----
    electron_mass_scaling = (base.m_p / base.m_e) / ratios_and_scaling.electron_proton_ratio
        #-------------------------------------------------------------------------------------
    charge_scaling = ratios_and_scaling.q_charge_scaling
    eps_0_scaling = ratios_and_scaling.eps_0_scaling
    mu_0_scaling = 1e0    
    scale_base_units(base, eps_0_scaling, mu_0_scaling, electron_mass_scaling, charge_scaling)
    #-----------------------------------------------


    #------------ physical constants depending on scaling ----------------
    mu = ( (base.m_e + base.m_p)  / 2. ) /  base.m_u# assuming only protons and electrons
    e = base.e # may be scaled
    c_real = base.c
    #------------------------------------------------------------------------
    
    k_B = base.k_B
    
    #--------- calculate real mass density ---------------
    number_density = known_real_units.number_density
    rho_real = number_density * 2 * mu * base.m_u
    #------------------------------------------------
    
    
    
    #------------ Load code units ----------------
    p_code = code_units.p
    rho_code = code_units.rho
    B_code = code_units.B
    c_code = code_units.c
    debye_code = code_units.debye
    #------------------------------------------------
    

    #------------- Calculate real values and scaling factors ----------------
    T = (p_code / rho_code ) * ( (mu * base.m_u) / k_B) * (c_real / c_code)^2
    p_real = ( rho_real / (mu * base.m_u) )  * k_B * T 
    mass_density_scale = rho_real / rho_code
    B_real = B_code * sqrt(p_real / p_code)
    length_scale = 1. / (debye_code * e) * sqrt(k_B * T / (4 * pi * number_density ))
    time_scale = length_scale * c_c / c_r 
    #----------------------------------------------------------------
    

    #------------- assume equal electron and proton temperature, and no temperature scale ----
    T_e = T
    T_p = T
    temperature_scale = 1e0 # don't scale temperature for now - can be added later
    #-----------------------------------------------------------------------------------------

    #---------- Assume commen length is the same as length scaling ---  
    length = length_scale
    #----------------------------------------------------------------


    #---------- create a new scaling object containing all the scaling factors ----
    scaling = ScalingCGS(base,
                number_density, length, T, T_e, T_e, B_real,
                length_scale, mass_density_scale, time_scale, temperature_scale)
    #--------------------------------------------------------------------------------
    return scaling
end




#c_u = code_units(c_c, B_c, p_c, rho_c, debye_c)
#r_u = known_real_units(1e9)
#eps_0_scaling = 1e2
#r_s = ratios_and_scaling(e_p_ratio, eps_0_scaling, charge_scaling)
#
#new_scaling = find_real_units(base, c_u, r_u, r_s)