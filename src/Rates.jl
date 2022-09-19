module Rates
	
import PhysicalConstants.CODATA2018: e, ε_0, ħ, k_B, m_e
using ProgressMeter;
using Unitful

include("Constants.jl");
include("ENZ.jl")

function normalisationB(q, d, m)
	n = 2 * m - 1
	ζ = n * π / d;
	return sqrt(
		ħ * ω_L / (ε_0 * ε_ρ * d) / (ζ^2 + q^2)
	)
end

function overlapΞ(q, Q, d, m)
	n = 2 * m - 1
	ζ = n * π / d;
	dQ = q - Q
	return dQ * (1 + exp(im * d * dQ)) / (dQ^2 - ζ^2)
end

function κ(q, Q, d, n_max, β_j)
	result = 0u"m^2/s"
	for i in 1:n_max
		B = normalisationB(q, d, i)
		Ξ = overlapΞ(q, Q, d, i)
		result += B * Ξ * β_j * e / ħ
	end
	result
end

function electronEnergy(k, k_z)
	ħ^2 * (k^2 + k_z^2) / (2. * m)
end

function electronDistribution(k, k_z, electronic_temperature)
	ε = electronEnergy(k, k_z);
	thermal = k_B * electronic_temperature;
	conduction_band_dos = (ħ^2 / (2. * π * m * k_B * electronic_temperature))^1.5;
	conduction_band_dos * exp( - ε / thermal) * 1u"1/m^3"
end

function transferredOutOfPlane(q, k, k_z, θ, ω)
	sqrt(k_z^2 + 2 * k * q * cos(θ) - q^2 - 2 * m * ω / ħ + 0.0u"1/m^2" * im)
end

function dk_dω(q, k, k_z, θ, ω_j)
	m / ħ / transferredOutOfPlane(q, k, k_z, θ, ω_j)
end

function gammaKernel(q, ω, k, k_z, dk_z, ω_j, β_j, d, n_max, electronic_temperature)
		2 * π * abs(κ(q, dk_z, d, n_max, β_j))^2 * electronDistribution(k, k_z, electronic_temperature) * ENZ.lorentzian_density_of_states(ω, ω_j, γ)
end

using Cubature;

function gammaKernelStripped(q, ω, k, k_z, θ, ω_j, β_j, d, n_max, electronic_temperature)
	k = k * 1u"1/m"
	k_z = k_z * 1u"1/m"
	ω = ω * 1u"1/s";
	dk_z = transferredOutOfPlane(q, k, k_z, θ, ω)
	if real(dk_z) == 0
		return 0
	else
		return	ustrip(gammaKernel(q, ω, k, k_z, dk_z, ω_j, β_j, d, n_max, electronic_temperature))
	end
end

function gammaIntegratedk(q, ω, ω_ENZ, d, electronic_temperature, n_max, p)
	rate = 0;
	error = 0;

	frequencies, hopfield = ENZ.polariton_eigenvalues_from_enz_ω(q, d, ω_ENZ, n_max);

	for j in 1:n_max+1
		ω_j = frequencies[j]
		β_j = hopfield[1, j]
	 	result = hcubature(
				x -> 2 * π * x[1] * gammaKernelStripped(q, ustrip(ω), x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0], [5e9, 2 * π, 5e9];
				reltol=1e-3
			)
		rate += result[1]
		error += result[2]
	end
	next!(p)
	rate, error
end


function gammaIntegratedkjω(q, ω, ω_ENZ, d, electronic_temperature, n_max, j, p)
	rate = 0;
	error = 0;

	frequencies, hopfield = ENZ.polariton_eigenvalues_from_enz_ω(q, d, ω_ENZ, n_max);

	
	ω_j = frequencies[j]
	β_j = hopfield[1, j]
 	result = hcubature(
				x -> 2 * π * x[1] * gammaKernelStrippedω(q, ustrip(ω), x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0], [5e9, 2 * π, 5e9];
				reltol=1e-3
		)
	rate += result[1]
	error += result[2]
	
	next!(p)
	rate, error
end

function gammaKernelStrippedω(q, ω, k, k_z, θ, ω_j, β_j, d, n_max, electronic_temperature)
	k = k * 1u"1/m"
	k_z = k_z * 1u"1/m"
	ω = ω * 1u"1/s";
	dk_z = transferredOutOfPlane(q, k, k_z, θ, ω)
	if real(dk_z) == 0
		return 0
	else
		dk_z = real(dk_z) + im * 1e-99u"1/m"
		return	real(ustrip(gammaKernel(q, ω, k, k_z, dk_z, ω_j, β_j, d, n_max, electronic_temperature) * m / ħ / dk_z))
	end
end

function gammaIntegratedωj(q, ω_ENZ, d, electronic_temperature, n_max, j, p)
	# println("Solving for", q, ω_ENZ)
	rate = 0;
	error = 0;

	frequencies, hopfield = ENZ.polariton_eigenvalues_from_enz_ω(q, d, ω_ENZ, n_max);

	ω_min = ustrip(ω_T);
	ω_max = ustrip(ω_L) * 1.1;

	ω_j = frequencies[j]
	β_j = hopfield[1, j]
	# result = hcubature(
	# 		x -> 2 * π * x[1] * gammaKernelStrippedω(q, x[4], x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0, ω_min], [1e9, 2 * π, 1e9, ω_max];
	# 		reltol=1e-3, maxevals=5000000
	# 	)
	result = hcubature(
			x -> 2 * π * x[1] * gammaKernelStrippedω(q, x[3], x[1], x[2], 0*π/2, ω_j, β_j, d, n_max, electronic_temperature), [0, 0, ω_min], [1e9, 1e9, ω_max];
			reltol=1e-3, maxevals=5000000
		)
	rate += result[1]
	error += result[2]
	next!(p)

	# println(rate, error)
	rate, error
end

function gammaIntegratedω(q, ω_ENZ, d, electronic_temperature, n_max, p)
	rate = 0;
	error = 0;

	frequencies, hopfield = ENZ.polariton_eigenvalues_from_enz_ω(q, d, ω_ENZ, n_max);

	ω_min = ustrip(ω_T);
	ω_max = 1.1 * ustrip(ω_L);

	for j in 1:n_max+1
		ω_j = frequencies[j]
		β_j = hopfield[1, j]
	 	result = hcubature(
				x -> 2 * π * x[1] * gammaKernelStrippedω(q, x[4], x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0, ω_min], [5e9, 2 * π, 5e9, ω_max];
				reltol=1e-2
			)
		rate += result[1]
		error += result[2]
	end
	next!(p)
	rate, error
end

function gammaKernelStrippedkz(q, k_zd, k, k_z, θ, ω_j, β_j, d, n_max, electronic_temperature)
	k = k * 1u"1/m"
	k_z = k_z * 1u"1/m"
	k_zd = k_zd * 1u"1/m";
	dk_z = abs(k_z - k_zd)
	ω = ħ / 2 / m * (k^2 + 2 * k * q * cos(θ) - q^2 - k_zd^2);

	if ω < ω_T
		return 0
	else
		return	real(ustrip(gammaKernel(q, ω, k, k_z, dk_z, ω_j, β_j, d, n_max, electronic_temperature)))
	end
end


function gammaIntegratedkzj(q, ω_ENZ, d, electronic_temperature, n_max, j, p)
	rate = 0;
	error = 0;

	frequencies, hopfield = ENZ.polariton_eigenvalues_from_enz_ω(q, d, ω_ENZ, n_max);

	ω_min = ustrip(ω_T);
	ω_max = 1.1 * ustrip(ω_L);

	ω_j = frequencies[j]
	β_j = hopfield[1, j]
 	result = hcubature(
			x -> 2 * π * x[1] * gammaKernelStrippedkz(q, x[4], x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0, 0], [1e9, 2 * π, 1e9, 1e9];
			reltol=1e-3
			#, maxevals=5000000
		)
	rate += result[1]
	error += result[2]
	next!(p)
	rate, error
end


function η_1(q, ω_ENZ, d, electronic_temperature, n_max, p)
	rate = 0;
	error = 0;

	frequencies, hopfield = ENZ.polariton_eigenvalues_from_enz_ω(q, d, ω_ENZ, n_max);

	ω_min = ustrip(ω_T);
	ω_max = 1.1 * ustrip(ω_L);

	results = zeros(n_max + 1)

	for j in 1:n_max+1
		ω_j = frequencies[j]
		β_j = hopfield[1, j]
	 	result = hcubature(
				x -> 2 * π * x[1] * gammaKernelStrippedω(q, x[4], x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0, ω_min], [1e9, 2 * π, 1e9, ω_max];
				reltol=1e-1
			)
		results[j] += result[1]*1u"1/s" / γ / ENZ.bose_einstein(ω_j, 300u"K");
	end
	next!(p)
	results
end

end