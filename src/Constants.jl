import PhysicalConstants.CODATA2018: e, c_0, ε_0, ħ, k_B, m_e
using Unitful

const β = 1u"1";
const β_L = 15.3e3u"m/s"
const ω_L = 1.8327e14u"1/s";
const ω_T = 0.818 * ω_L;
const γ = ω_L / 240.;
const ε_∞ = 6.52;
const ε_st = 9.7;
const ε_ρ = 1. / (1. / ε_∞ - 1. / ε_st);
const ε_c = 11.71;
const m = 0.29 * m_e;
const k_min = 1e6u"1/m";
const k_max = 1e8u"1/m";
