import PhysicalConstants.CODATA2018: m_e

const β_L::Float64 = 15.3e3;
const ω_L::Float64 = 1.8327e14;
const ω_T::Float64 = 0.818 * ω_L;
const γ::Float64 = ω_L / 240.;
const ε_∞::Float64 = 6.52;
const ε_st::Float64 = 9.7;
const ε_ρ::Float64 = 1. / (1. / ε_∞ - 1. / ε_st);
const ε_c::Float64 = 11.71;
const m::Float64 = 0.29 * ustrip(m_e);
const e::Float64 = 1.9e-19;
const ħ::Float64 = 1.064e-34;
const k_B::Float64 = 1.380649e-23;
const ε_0::Float64 = 8.854e-12;