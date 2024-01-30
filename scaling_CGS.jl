include("Base_units.jl")

"""
    ScalingCGS

A struct for handling scaling in the CGS (Centimeter-Gram-Second) unit system.
"""
mutable struct ScalingCGS
    system::String  # Unit system ('CGS' in this case)

    # Input physical parameters
    number_density_real::Float64       # Number density [cm^-3]
    length_real::Float64               # Characteristic length [cm]
    temperature_real::Float64          # Temperature [K]
    B_flux_real::Float64               # Magnetic flux density [Gauss]

    # Scaling factors
    length_scaling::Float64            # Length scaling factor
    mass_density_scaling::Float64      # Mass density scaling factor
    time_scaling::Float64              # Time scaling factor

    # Base physical constants (from BaseUnits)
    base_units::BaseUnits

    # Derived real values
    mu_real::Float64                   # Mean molecular weight [unitless]
    mass_density_real::Float64         # Mass density [g/cm^3]
    time_real::Float64                 # Time [s]
    u_real::Float64                    # Characteristic velocity [cm/s]
    c_real::Float64                    # Speed of light [cm/s]
    v_a_real::Float64                  # Alfven velocity [cm/s]
    v_thermal_e_real::Float64          # Electron thermal velocity [cm/s]
    v_thermal_p_real::Float64          # Proton thermal velocity [cm/s]
    electron_gyro_freq_real::Float64   # Electron gyro frequency [rad/s]
    proton_gyro_freq_real::Float64     # Proton gyro frequency [rad/s]
    electron_plasma_freq_real::Float64 # Electron plasma frequency [rad/s]
    proton_plasma_freq_real::Float64   # Proton plasma frequency [rad/s]
    skin_depth_real::Float64           # Electron skin depth [cm]
    skin_depth_real_p::Float64         # Proton skin depth [cm]
    debye_len_real::Float64            # Debye length [cm]
    inter_dist_real::Float64           # Interparticle distance [cm]
    e_gyro_radiues_real::Float64       # Electron gyro radius [cm]
    p_gyro_radiues_real::Float64       # Proton gyro radius [cm]
    mass_real::Float64                 # Characteristic mass [g]
    energy_real::Float64               # Characteristic energy [erg]
    energy_per_mass_real::Float64      # Energy per mass [erg/g]
    energy_per_volume_real::Float64    # Energy per volume [erg/cm^3]
    pressure_real::Float64             # Pressure [Ba]

    # Additional fields for scaling factors and other properties can be added here

    #scaling values 
    mu_scaling :: Float64
    number_density_scaling::Float64
    u_scaling::Float64
    c_scaling::Float64
    v_a_scaling::Float64
    v_thermal_e_scaling::Float64
    v_thermal_p_scaling::Float64
    electron_gyro_freq_scaling::Float64
    proton_gyro_freq_scaling::Float64
    electron_plasma_freq_scaling::Float64
    proton_plasma_freq_scaling::Float64
    skin_depth_scaling::Float64
    skin_depth_p_scaling::Float64
    debye_len_scaling::Float64
    inter_dist_scaling::Float64
    e_gyro_radiues_scaling::Float64
    p_gyro_radiues_scaling::Float64
    mass_scaling::Float64
    B_flux_scaling :: Float64
    temperature_scaling::Float64
    energy_scaling::Float64
    energy_per_mass_scaling::Float64
    energy_per_volume_scaling::Float64
    pressure_scaling::Float64

    #code values 
    mu_code :: Float64
    mass_density_code::Float64
    number_density_code::Float64
    u_code::Float64
    c_code::Float64
    v_a_code::Float64
    v_thermal_e_code::Float64
    v_thermal_p_code::Float64
    time_code :: Float64
    electron_gyro_freq_code::Float64
    proton_gyro_freq_code::Float64
    electron_plasma_freq_code::Float64
    proton_plasma_freq_code::Float64
    length_code :: Float64
    skin_depth_code::Float64
    skin_depth_p_code::Float64
    debye_len_code::Float64
    inter_dist_code::Float64
    e_gyro_radiues_code::Float64
    p_gyro_radiues_code::Float64
    mass_code::Float64
    B_flux_code :: Float64
    temperature_code::Float64
    energy_code::Float64
    energy_per_mass_code::Float64
    energy_per_volume_code::Float64
    pressure_code::Float64

    function ScalingCGS(base_units, number_density=1e15, length=1e8, temperature=1e6, B_flux=1e1, length_scale=3e8, mass_density_scale=1e-13, time_scale=1.0)
        system = base_units.system
        obj = new(system, number_density, length, temperature, B_flux, length_scale, mass_density_scale, time_scale, base_units)
        # Additional initialization and calculations for derived values
        return obj
    end
end




# Example usage
base = BaseUnits("CGS")
scaling = ScalingCGS(base)
