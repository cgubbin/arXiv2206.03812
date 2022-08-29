include("../src/Rates.jl")
include("../src/ENZ.jl")

using .ENZ;
using .Rates;
using CairoMakie;
using Interpolations;
using Unitful;
using ThreadsX;
using ProgressMeter;


# # instantiate and precompile environment in all processes
# @everywhere begin
#   using Pkg; Pkg.activate(@__DIR__)
#   Pkg.instantiate(); Pkg.precompile()
# end

import PhysicalConstants.CODATA2018: c_0;

function wavenumber_to_omega(wn)
	k = wn * 2. * π;
	return k * c_0
end


d = 2e-9u"m"


# Precompute the epsilon_near_zero_dispersion
minimum_wavevector = ENZ.ω_L * sqrt(ENZ.ε_c) / ENZ.c_0 * 1.01;
maximum_wavevector = 1e9u"1/m"
wavevectors = LinRange(minimum_wavevector, maximum_wavevector, 5000);
initial = [1.75];
enz_frequencies = [ENZ.epsilon_near_zero_dispersion!(k, d, initial) for (idx, k) in enumerate(wavevectors)]

ω_ENZ_int = linear_interpolation(wavevectors, enz_frequencies)

q_x = 1e8u"1/m"

q_min = minimum_wavevector;
q_max = 1e8u"1/m"
n_bin_k = 200;

n_bin_ω = 200;
wavevectors = LinRange(minimum_wavevector, maximum_wavevector, n_bin_k);

# wavenumbers = LinRange(82500u"1/m", 100000u"1/m", n_bin_ω);
# frequencies = [wavenumber_to_omega(wn) for wn in wavenumbers]

	
# n = n_bin_k * n_bin_ω
# p = Progress(n, 1)

# map = ThreadsX.collect(Rates.gammaIntegratedk(
# 		q, ω, ω_ENZ_int(q), 2e-9u"m", 1, 300u"K", 2, p
# 	)[1] * 1e53 for q in wavevectors, ω in frequencies
# )

# fig = Figure()
# ax, hm = heatmap(
# 	fig[1, 1][1, 1], ustrip(wavevectors), ustrip(wavenumbers), map,
# 	)
# Colorbar(fig[1, 1][1, 2], hm)
# fig
# map

n = n_bin_k
p = Progress(n, 1)

# results = ThreadsX.collect(Rates.gammaIntegratedω(
# 		q, ω_ENZ_int(q), 2e-9u"m", 1, 300u"K", 2, p
# 	)[1] * 1e53 for q in wavevectors
# )

# results

q = wavevectors[10]
T = 700u"K"

# Rates.gammaIntegratedω(
# 		q, ω_ENZ_int(q), 2e-9u"m", T, 2, p
# 	)
	
	#[1] / ENZ.γ * 1e24 / ENZ.bose_einstein(ω_ENZ_int(q), 750u"K")
	


Rates.η_1(
		q, ω_ENZ_int(q), 2e-9u"m", T, 2, p
	)