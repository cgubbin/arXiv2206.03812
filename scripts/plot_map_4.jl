
include("../src/ENZ.jl")

using CairoMakie;
using DelimitedFiles
using .ENZ;
using Interpolations;
# using LaTeXStrings;
using TensorCast;
using Unitful;
using FStrings;
using ColorSchemes;

import PhysicalConstants.CODATA2018: c_0;

function wavenumber_to_omega(wn)
	k = wn * 2. * π;
	return k * c_0
end

thickness_target = "2nm";
d = 2e-9u"m";
n_max = 2;
files_to_act_on = filter(contains("2nm"), readdir("results"))

### Act to construct the epsilon_near_zero_dispersion
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
@cast pol_frequencies[i,j] := pol_frequencies[i][j]
pol_eigenvectors = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[2] for k in wavevectors]
@cast pol_eigenvectors[i,j,k] := pol_eigenvectors[i][j,k]

T_lattice = 300u"K"

n_thermal = [ENZ.bose_einstein(ω, T_lattice) for ω in pol_frequencies]

for file in files_to_act_on
	println(file)
	mat = readdlm("results/" * file, ',', Float64);

	electrical_map = [
		sum([mat[j, m] * ENZ.lorentzian_density_of_states(ω, pol_frequencies[j, m], ENZ.γ) for m in 1:n_max+1])
		for j in 1:n_bin_k, ω in frequencies 
	]

	electrical_map_vg = [
		sum([abs(vgs(wavevectors[j])[m][1]) * mat[j, m] * ENZ.lorentzian_density_of_states(ω, pol_frequencies[j, m], ENZ.γ) for m in 1:n_max+1])
		for j in 1:n_bin_k, ω in frequencies 
	]

	thermal_map = [
		sum([ENZ.bose_einstein(pol_frequencies[j, m], T_lattice) * ENZ.lorentzian_density_of_states(ω, pol_frequencies[j, m], ENZ.γ) for m in 1:n_max+1])
		for j in 1:n_bin_k, ω in frequencies 
	]

	thermal_map_vg = [
		sum([abs(vgs(wavevectors[j])[m][1]) * ENZ.bose_einstein(pol_frequencies[j, m], T_lattice) * ENZ.lorentzian_density_of_states(ω, pol_frequencies[j, m], ENZ.γ) for m in 1:n_max+1])
		for j in 1:n_bin_k, ω in frequencies 
	]




	f = Figure(resolution = (1400, 1100), font = "CMU Serif", xlabelsize=22)

	ax = heatmap(f[1, 1][1, 1], ustrip(wavevectors) / 1e9, ustrip(wavenumbers) / 100, ustrip(thermal_map) / 1e-13, colorrange=(0, 1),
		axis = (title = L"\mathrm{Thermal} \; \times \; 10", xlabel = L"\mathrm{LTP\;Wavevector}\;(\mathrm{nm}^{-1})", ylabel = L"\mathrm{Wavenumber}\;(\mathrm{cm}^{-1})"), interpololate = true,
		colormap = :lajolla
	)
	heatmap(f[1, 1][1, 2], ustrip(wavevectors) / 1e9, ustrip(wavenumbers) / 100, ustrip(thermal_map_vg)  / 1e-17, colorrange=(0, 1),
		axis = (title = L"\mathrm{Thermal}\;\times \; 10 \; \nu_g / c", xlabel = L"\mathrm{LTP\;Wavevector}\;(\mathrm{nm}^{-1})", ylabel =  L"\mathrm{Wavenumber}\;(\mathrm{cm}^{-1})"), interpololate = true,
		colormap = :lajolla
	)

	heatmap(f[1, 1][2, 1], ustrip(wavevectors) / 1e9, ustrip(wavenumbers) / 100, ustrip(electrical_map) / 1e-13, colorrange=(0, 10),
		axis = (title = L"\mathrm{Electrical}", xlabel = L"\mathrm{LTP\;Wavevector}\;(\mathrm{nm}^{-1})", ylabel = L"\mathrm{Wavenumber}\;(\mathrm{cm}^{-1})"), interpololate = true,
		colormap = :lajolla
	)
	heatmap(f[1, 1][2, 2], ustrip(wavevectors) / 1e9, ustrip(wavenumbers) / 100, ustrip(electrical_map_vg)   / 1e-17, colorrange=(0, 10),
		axis = (title = L"\mathrm{Electrical}\; \times \; \nu_g / c", xlabel = L"\mathrm{LTP\;Wavevector}\;(\mathrm{nm}^{-1})", ylabel =  L"\mathrm{Wavenumber}\;(\mathrm{cm}^{-1})"), interpololate = true,
		colormap = :lajolla
	)

# hidedecorations!.([ax, current_axis()])
# ax.ylabelvisible = true

# f[1, 1][2, 1:2, Bottom()] = Label(f, L"\mathrm{LTP\;Wavevector}\;(\mathrm{nm}^{-1})", padding = (0, 0, 0, 10))
f[1, 1][3, 1] = Colorbar(f, height = 20, vertical = false, label = "", colormap = :lajolla,
  # ticklabelalign = (:center, :top), 
	limits = (0, 1e-13),
	 ticks=[0.0, 0.5e-13, 1e-13])
f[1, 1][3, 2] = Colorbar(f, height = 20, vertical = false, label = "", colormap = :lajolla,
  # ticklabelalign = (:center, :top),
	limits = (0, 1e-17),
	 ticks=[0.0, 0.5e-17, 1e-17])

	substrings = split(split(file, ".")[1], "_")
	t_term = substrings[3]
	d_term = substrings[2]

	save(f"figures/maps_{d_term}_{t_term}.pdf", current_figure())
end;

