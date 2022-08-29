include("../src/ENZ.jl")

using .ENZ;
using CairoMakie;
using Interpolations;
using Unitful;
using LaTeXStrings




d = 1.5e-9u"m"

minimum_wavevector = ENZ.ω_L * sqrt(ENZ.ε_c) / ENZ.c_0 * 1.01;
maximum_wavevector = 1e9u"1/m"
wavevectors = LinRange(minimum_wavevector, maximum_wavevector, 1000);


initial = [1.75];
enz_frequencies = [ENZ.epsilon_near_zero_dispersion!(k, d, initial) for (idx, k) in enumerate(wavevectors)]

enz_interp = linear_interpolation(wavevectors, enz_frequencies)
interp = [enz_interp(k) for k in wavevectors];# * 1u"1/s";

initial = [1.75];
pol_frequencies = [ENZ.polariton_eigenvalues!(k, d, initial, 2)[1] for k in wavevectors]
pol_frequencies = mapreduce(permutedims, vcat, pol_frequencies);

minimum_frequency = ENZ.ω_T
maximum_frequency = ENZ.ω_L * 1.1
frequencies = LinRange(minimum_frequency, maximum_frequency, 1000)

initial = [1.75]
density_of_states = [ENZ.polariton_density_of_states(ω, k, d, initial, 2) for k in wavevectors, ω in frequencies]
density_of_states /= maximum(density_of_states)

fig, ax, hm = heatmap(
	ustrip(wavevectors), ustrip(ENZ.omega_to_wavenumber(frequencies)), ustrip(density_of_states),
	colormap = :lajolla, colorrange = (0, 1)
)
ax.xlabel = "Parallel Wavevector 1/m"
ax.ylabel = "Wavenumber 1/m"
Colorbar(fig[:, end+1], hm)
lines!(ustrip(wavevectors), ustrip(ENZ.omega_to_wavenumber(interp)), color = :green, linewidth=4)
lines!(ustrip(wavevectors), ustrip(ENZ.omega_to_wavenumber(pol_frequencies[:, 1])), color = :red, linewidth=4)
lines!(ustrip(wavevectors), ustrip(ENZ.omega_to_wavenumber(pol_frequencies[:, 2])), color = :red, linewidth=4)
lines!(ustrip(wavevectors), ustrip(ENZ.omega_to_wavenumber(pol_frequencies[:, 3])), color = :red, linewidth=4)
ylims!(low=82500, high=100000)
xlims!(low=0, high=ustrip(maximum_wavevector))
save("2nm.pdf", fig)

fig
