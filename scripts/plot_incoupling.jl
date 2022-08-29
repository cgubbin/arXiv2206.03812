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


# Precompute the epsilon_near_zero_dispersion after checking in plot_dispersion.jl
minimum_wavevector = ENZ.ω_L * sqrt(ENZ.ε_c) / ENZ.c_0 * 1.01;
maximum_wavevector = 1e9u"1/m"
wavevectors = LinRange(minimum_wavevector, maximum_wavevector, 5000);
initial = [1.75];
enz_frequencies = [ENZ.epsilon_near_zero_dispersion!(k, d, initial) for (idx, k) in enumerate(wavevectors)]

## The interpolated ENZ dispersion over the range of interest
ω_ENZ_int = linear_interpolation(wavevectors, enz_frequencies)

## Generate the plotting grid
q_min = minimum_wavevector;
q_max = 1e8u"1/m"
n_bin_k = 200;
n_bin_ω = 200;
wavevectors = LinRange(minimum_wavevector, maximum_wavevector, n_bin_k);
wavenumbers = LinRange(82500u"1/m", 100000u"1/m", n_bin_ω);
frequencies = [wavenumber_to_omega(wn) for wn in wavenumbers]

	
n = n_bin_k * n_bin_ω
p = Progress(n, 1)

T = 800u"K"

map = ThreadsX.collect(Rates.gammaIntegratedk(
		q, ω, ω_ENZ_int(q), 2e-9u"m", T, 2, p
	)[1] * 1u"s" / ENZ.γ * 1e24 for q in wavevectors, ω in frequencies
)

peak = maximum(map);
map /= peak;

fig = Figure()
ax, hm = heatmap(
	fig[1, 1][1, 1], ustrip(wavevectors), ustrip(wavenumbers), map,
	colormap = :lajolla, colorrange = (0, 1), interpolate=true
	)
Colorbar(fig[1, 1][1, 2], hm)
ax.xlabel = "Parallel Wavevector 1/m"
ax.ylabel = "Wavenumber 1/m"
save("2nm_800K.pdf", fig)
fig