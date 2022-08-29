include("../src/Rates.jl")
include("../src/ENZ.jl")

using .ENZ;
using .Rates;
using CairoMakie;
using Interpolations;
using Unitful;
using ThreadsX;
using ProgressMeter;

import PhysicalConstants.CODATA2018: c_0;

function wavenumber_to_omega(wn)
	k = wn * 2. * π;
	return k * c_0
end


d = 2e-9u"m"
n_max = 2;
T_lattice = 300u"K"


# Precompute the epsilon_near_zero_dispersion
minimum_wavevector = ENZ.ω_L * sqrt(ENZ.ε_c) / ENZ.c_0 * 1.4;
maximum_wavevector = 1e9u"1/m"
wavevectors = LinRange(minimum_wavevector, maximum_wavevector, 5000);
initial = [1.75];
enz_frequencies = [ENZ.epsilon_near_zero_dispersion!(k, d, initial) for (idx, k) in enumerate(wavevectors)]

ω_ENZ_int = linear_interpolation(wavevectors, enz_frequencies)

initial = [1.75];
pol_frequencies = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[1] for k in wavevectors]
pol_frequencies = mapreduce(permutedims, vcat, pol_frequencies);

ω_p1_int = linear_interpolation(wavevectors, pol_frequencies[:, 1])
ω_p2_int = linear_interpolation(wavevectors, pol_frequencies[:, 2])
ω_p3_int = linear_interpolation(wavevectors, pol_frequencies[:, 3])

vg1(x) = Interpolations.gradient(ω_p1_int, x) / c_0
vg2(x) = Interpolations.gradient(ω_p2_int, x) / c_0
vg3(x) = Interpolations.gradient(ω_p3_int, x) / c_0

vgs(x) = [vg1(x), vg2(x), vg3(x)]

n_bin_k = 200;
n_bin_ω = 200;
wavevectors = LinRange(minimum_wavevector, maximum_wavevector, n_bin_k);
wavenumbers = LinRange(82500u"1/m", 100000u"1/m", n_bin_ω);
frequencies = [wavenumber_to_omega(wn) for wn in wavenumbers]

initial = [1.75];
pol_frequencies = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[1] for k in wavevectors]
pol_frequencies = mapreduce(permutedims, vcat, pol_frequencies);

thermal_map = [
	sum([ENZ.bose_einstein(pol_frequencies[j, m], T_lattice) * ENZ.lorentzian_density_of_states(ω, pol_frequencies[j, m], ENZ.γ) for m in 1:n_max+1])
	for j in 1:n_bin_k, ω in frequencies 
]

thermal_map_vg = [
	sum([abs(vgs(wavevectors[j])[m][1]) * ENZ.bose_einstein(pol_frequencies[j, m], T_lattice) * ENZ.lorentzian_density_of_states(ω, pol_frequencies[j, m], ENZ.γ) for m in 1:n_max+1])
	for j in 1:n_bin_k, ω in frequencies 
]

fig = Figure(resolution=(1200, 500))
ax1, hm1 = heatmap(
	fig[1, 1][1, 1], ustrip(wavevectors / 100), ustrip(wavenumbers / 100), ustrip(thermal_map),
	colormap = :lajolla, interpolate=true, colorrange = (0, 1.001e-13)
	)
ax2, hm2 = heatmap(
	fig[1, 2][1, 1], ustrip(wavevectors / 100), ustrip(wavenumbers / 100), ustrip(thermal_map_vg),
	colormap = :lajolla, interpolate=true, colorrange = (0, 1e-17)
	)
Colorbar(fig[1, 1][1, 2], hm1)
Colorbar(fig[1, 2][1, 2], hm2)

ax1.xlabel = "Parallel Wavevector 1/cm"
ax2.xlabel = "Parallel Wavevector 1/cm"
ax1.ylabel = "Wavenumber 1/cm"
ax2.ylabel = "Wavenumber 1/cm"
save("figure_4_2nm_300K.pdf", fig)
fig