module FieldBalance

include("Constants.jl");

import PhysicalConstants.CODATA2018: e, ε_0, ħ, k_B, m_e
import SpecialFunctions: besselk
using Unitful

const E_0 = m * e * ħ * ω_L / (ħ^2 * ε_0) * (1. / ε_∞ - 1. / ε_st);

const Θ_D = ħ * ω_L / k_B;


function bose_einstein_distribution(ω, T)
	return 1 / (exp(ħ * ω / (k_B * T)) - 1.)
end

function drift_velocity(electronic_temperature, lattice_temperature)
 	x_o = ħ * ω_L / (k_B * lattice_temperature);
 	x_e = ħ * ω_L / (k_B * electronic_temperature);

	return sqrt(
	3. * k_B * Θ_D / (m * x_e) 
	* (exp(x_o - x_e) - 1.) * besselk(0, x_e / 2.) / ((exp(x_o - x_e) + 1) * besselk(1, x_e / 2.) + (exp(x_o - x_e) - 1) * besselk(0, x_e / 2.))
	)
end

function electric_field(electronic_temperature, lattice_temperature)
 	x_o = ħ * ω_L / (k_B * lattice_temperature);
 	x_e = ħ * ω_L / (k_B * electronic_temperature);
	n = bose_einstein_distribution(ω_L, lattice_temperature);

	return E_0 * sqrt(
		2. / (3. * π) * n^2 * x_e^2 * exp(x_e) * besselk(0, x_e / 2.) * (
		(exp(x_o - x_e) + 1) * besselk(1, x_e / 2.) + (exp(x_o - x_e) - 1) * besselk(0, x_e / 2.)
		)
	)
end

using Cubature;

function bulk_like_absorption_percentage(electronic_temperature, lattice_temperature, d)
	ν_D = drift_velocity(electronic_temperature, lattice_temperature);



	(res, err) = hcubature(x -> begin bulkIntegrand(x); end, [0,0],[1,1])
	Γ = res * 1u"1/s";
	α = Γ / ν_D;
	return (1 - exp(- α * d)) * 100.
end

function layer_like_absorption_percentage(electronic_temperature, lattice_temperature, d)
	ν_D = drift_velocity(electronic_temperature, lattice_temperature);



	(res, err) = hcubature(x -> begin layerIntegrand(x); end, [0,0,0],[1,1,1])
	Γ = res * 1u"1/s";
	α = Γ / ν_D;
	return (1 - exp(- α * d)) * 100.
end

# println(bulk_like_absorption_percentage(1000u"K", 300u"K", 1e-9u"m"))
# println(layer_like_absorption_percentage(1000u"K", 300u"K", 1e-9u"m"))

end