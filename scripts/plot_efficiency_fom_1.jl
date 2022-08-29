include("../src/ENZ.jl")

using CairoMakie;
using DelimitedFiles
using .ENZ;
using Interpolations;
# using LaTeXStrings;
using TensorCast;
using Unitful;
using FStrings;

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
	mat = readdlm("results/" * file, ',', Float64);
  lines(ustrip(wavevectors), mat[:, 1] ./ n_thermal[:, 1]; color = "#56B4E9", linewidth = 2, label = L"\mathrm{Branch}\;1",
        axis = (xlabel = L"\mathrm{LTP\;Wavevector}\;(\mathrm{m}^{-1})", ylabel = L"\mathrm{Figure\;of\;Merit\;\;}\eta^{(1)}", xgridcolor = :red,
            xlabelsize = 22, ylabelsize = 22,
            xgridstyle = :dashdot, xgridwidth = 0.85,
            xtickalign = 1, xticksize = 20),
        figure = (resolution = (600, 400), font = "CMU Serif"))
	# break
  lines!(ustrip(wavevectors), mat[:, 2] ./ n_thermal[:, 2]; color = :black, linestyle = :dash, label = L"\mathrm{Branch}\;2")
  lines!(ustrip(wavevectors), mat[:, 3] ./ n_thermal[:, 3]; color = :red, linestyle = :dash, label = L"\mathrm{Branch}\;3")
  xlims!(minimum(ustrip(wavevectors)), maximum(ustrip(wavevectors)))
	ylims!(0, 50)
  axislegend("Legend", position = :lt)
	substrings = split(split(file, ".")[1], "_")
	t_term = substrings[3]
	d_term = substrings[2]

	save(f"figures/efficiency_{d_term}_{t_term}.pdf", current_figure())
end;
