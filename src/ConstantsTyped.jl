import PhysicalConstants.CODATA2018: m_e

const β_L::Float64 = 15.3e3; # The velocity of an LO phonon in m/s
const ω_L::Float64 = 1.8327e14; # The longitudinal phonon frequency in rad/s
const ω_T::Float64 = 0.818 * ω_L; # The transverse phonon frequency in rad/s
const γ::Float64 = ω_L / 240.; # The phonon damping in rad / s
const ε_∞::Float64 = 6.52; # The high-frequency dielectric constant
const ε_st::Float64 = 9.7; # The static dielectric constant
const ε_ρ::Float64 = 1. / (1. / ε_∞ - 1. / ε_st); # The Fröhlich constant
const ε_c::Float64 = 11.71; # The cladding dielectric function (represents Si)
const m::Float64 = 0.29 * ustrip(m_e); # The effective mass in the SiC
const e::Float64 = 1.9e-19; # The unit charge
const ħ::Float64 = 1.064e-34; # Reduced Planck constant in J s
const k_B::Float64 = 1.380649e-23; # The Boltzmann constant in J / K
const ε_0::Float64 = 8.854e-12; # The permitivitty of free-space