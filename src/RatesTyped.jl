module RatesTyped
	
using ProgressMeter;
using Unitful

include("ConstantsTyped.jl");
include("ENZ.jl")

function normalisationB(q::Float64, d::Float64, m::Int64)
	n::Int64 = 2 * m - 1
	ζ::Float64 = n * π / d;
	return sqrt(
		ħ * ω_L / (ε_0 * ε_ρ * d) / (ζ^2 + q^2)
	)
end

function overlapΞ(q::Float64, Q::Float64, d::Float64, m::Int64)
	n::Int64 = 2 * m - 1
	ζ::Float64 = n * π / d;
	dQ::Float64 = q - Q
	return dQ * (1 + exp(im * d * dQ)) / (dQ^2 - ζ^2)
end

function κ(q::Float64, Q::Float64, d::Float64, n_max::Int64, β_j::ComplexF64)
	result::ComplexF64 = 0.;
	for i in 1:n_max
		B::Float64 = normalisationB(q, d, i)
		Ξ::ComplexF64 = overlapΞ(q, Q, d, i)
		result += B * Ξ * β_j * e / ħ
	end
	result
end

function electronEnergy(k::Float64, k_z::Float64)
	ħ^2 * (k^2 + k_z^2) / (2. * m)
end

function electronDistribution(k::Float64, k_z::Float64, electronic_temperature::Float64)
	ε::Float64 = electronEnergy(k, k_z);
	thermal::Float64 = k_B * electronic_temperature;
	conduction_band_dos::Float64 = (ħ^2 / (2. * π * m * k_B * electronic_temperature))^1.5;
	conduction_band_dos * exp( - ε / thermal)
end

function transferredOutOfPlane(q::Float64, k::Float64, k_z::Float64, θ::Float64, ω::Float64)
	sqrt(k_z^2 + 2 * k * q * cos(θ) - q^2 - 2 * m * ω / ħ + 0.0 * im)
end

function dk_dω(q::Float64, k::Float64, k_z::Float64, θ::Float64, ω_j::Float64)
	m / ħ / transferredOutOfPlane(q, k, k_z, θ, ω_j)
end

function gammaKernel(q::Float64, ω::Float64, k::Float64, k_z::Float64, dk_z::Float64, ω_j::Float64, β_j::ComplexF64, d::Float64, n_max::Int64, electronic_temperature::Float64)
		2 * π * abs(κ(q, dk_z, d, n_max, β_j))^2 * electronDistribution(k, k_z, electronic_temperature) * ENZ.lorentzian_density_of_states(ω, ω_j, γ)
end

function gammaKernelAbsorption(q::Float64, ω::Float64, k::Float64, k_z::Float64, dk_z::Float64, ω_j::Float64, β_j::ComplexF64, d::Float64, n_max::Int64, electronic_temperature::Float64)
		2 * π * abs(κ(q, dk_z, d, n_max, β_j))^2 * electronDistribution(k, k_z, electronic_temperature) * ENZ.lorentzian_density_of_states(ω, -ω_j, γ)
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

function gammaKernelStrippedω(q::Float64, ω::Float64, k::Float64, k_z::Float64, θ::Float64, ω_j::Float64, β_j::ComplexF64, d::Float64, n_max::Int64, electronic_temperature::Float64)
	# k = k * 1u"1/m"
	# k_z = k_z * 1u"1/m"
	# ω = ω * 1u"1/s";
	dk_z::ComplexF64 = transferredOutOfPlane(q, k, k_z, θ, ω)
	if real(dk_z) == 0
		return 0
	else
		dk_z = real(dk_z) + im * 1e-99
		return	real(gammaKernel(q, ω, k, k_z, real(dk_z), ω_j, β_j, d, n_max, electronic_temperature) * m / ħ / dk_z)
	end
end

function gammaIntegratedωj(q::Float64, ω_ENZ::Float64, d::Float64, electronic_temperature::Float64, n_max::Int64, j::Int64, p::Progress)
	# println("Solving for", q, ω_ENZ)
	rate::Float64 = 0;
	error::Float64 = 0;

	frequencies, hopfield = ENZ.polariton_eigenvalues_from_enz_ω(q*1u"1/m", d*1u"m", ω_ENZ*1u"1/s", n_max);

	ω_min::Float64 = 1.0 * ω_T;
	ω_max::Float64 = 1.1 * ω_L;

	ω_j::Float64 = ustrip(frequencies[j])
	β_j::ComplexF64 = hopfield[1, j]
	result = pcubature(
			x -> 2 * π * x[1] * gammaKernelStrippedω(q, x[4], x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0, ω_min], [5e9, 2 * π, 5e9, ω_max];
			reltol=1e-3
		)
	rate += result[1]
	error += result[2]
	next!(p)

	# println(rate, error)
	rate, error
end

function gammaKernelStrippedkz(q::Float64, k_zd::Float64, k::Float64, k_z::Float64, θ::Float64, ω_j::Float64, β_j::ComplexF64, d::Float64, n_max::Int64, electronic_temperature::Float64)
	dk_z::Float64 = abs(k_z - k_zd)
	ω::Float64 = ħ / 2 / m * (k_z^2 + 2 * k * q * cos(θ) - q^2 - k_zd^2);

	if ω < ω_T || ω > 1.1 * ω_L
		return 0
	else
		return	real(ustrip(gammaKernel(q, ω, k, k_z, dk_z, ω_j, β_j, d, n_max, electronic_temperature)))
	end
end

function gammaKernelStrippedAbsorptionkz(q::Float64, k_zd::Float64, k::Float64, k_z::Float64, θ::Float64, ω_j::Float64, β_j::ComplexF64, d::Float64, n_max::Int64, electronic_temperature::Float64)
	dk_z::Float64 = abs(k_z - k_zd)
	ω::Float64 = ħ / 2 / m * (k_z^2 - 2 * k * q * cos(θ) - q^2 - k_zd^2);

	if ω > -ω_T || ω < - 1.1 * ω_L
		return 0
	else
		return	real(ustrip(gammaKernelAbsorption(q, ω, k, k_z, dk_z, ω_j, β_j, d, n_max, electronic_temperature)))
	end
end

function gammaIntegratedAbsorptionkzj(q, ω_pol, β_enz, d, electronic_temperature, n_max, j)

	rate::Float64 = 0;
	error::Float64 = 0;


	ω_j::Float64 = ustrip(ω_pol[j](q));
	β_j::ComplexF64 = β_enz[j](q);

	q = ustrip(q);
	electronic_temperature = ustrip(electronic_temperature);
	d = ustrip(d);
	 
 	result = hcubature(
			x -> 2 * π * x[1] * gammaKernelStrippedAbsorptionkz(q, x[4], x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0, 0], [1e9, 2 * π, 1e9, 1e9];
			reltol=1e-2, maxevals=5000000
		)
	rate += result[1]
	error += result[2]
	# next!(p)
    # If we timed out then return 0
    if abs(error) / abs(rate) > 1e-2
        rate = 0.0
    end
	rate, error
end

function gammaIntegratedkzj(q, ω_pol, β_enz, d, electronic_temperature, n_max, j)

	rate::Float64 = 0;
	error::Float64 = 0;


	ω_j::Float64 = ustrip(ω_pol[j](q));
	β_j::ComplexF64 = β_enz[j](q);

	q = ustrip(q);
	electronic_temperature = ustrip(electronic_temperature);
	d = ustrip(d);
	 
 	result = hcubature(
			x -> 2 * π * x[1] * gammaKernelStrippedkz(q, x[4], x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0, 0], [1e9, 2 * π, 1e9, 1e9];
			reltol=1e-2, maxevals=5000000
		)
	rate += result[1]
	error += result[2]
	# next!(p)
    # If we timed out then return 0
    if abs(error) / abs(rate) > 1e-2
        rate = 0.0
    end
	rate, error
end

function gammaIntegratedkzjo(q, ω_ENZ, d, electronic_temperature, n_max, j, p)

	rate::Float64 = 0;
	error::Float64 = 0;

	frequencies, hopfield = ENZ.polariton_eigenvalues_from_enz_ω(q*1u"1/m", d*1u"m", ω_ENZ*1u"1/s", n_max);

	ω_min::Float64 = 1.0 * ω_T;
	ω_max::Float64 = 1.1 * ω_L;

	ω_j::Float64 = ustrip(frequencies[j])
	β_j::ComplexF64 = hopfield[1, j]
	 
 	result = hcubature(
			x -> 2 * π * x[1] * gammaKernelStrippedkz(q, x[4], x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0, 0], [1e9, 2 * π, 1e9, 1e9];
			reltol=1e-2# maxevals=5000000
		)
	rate += result[1]
	error += result[2]
	next!(p)
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
				x -> 2 * π * x[1] * gammaKernelStrippedω(q, x[4], x[1], x[3], x[2], ω_j, β_j, d, n_max, electronic_temperature), [0, 0, 0, ω_min], [5e9, 2 * π, 5e9, ω_max];
				reltol=1e-2
			)
		results[j] += result[1]*1u"1/s" / γ / ENZ.bose_einstein(ω_j, 300u"K");
	end
	next!(p)
	results
end

end
