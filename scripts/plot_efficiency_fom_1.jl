include("../src/ENZ.jl")
include("EvaluationConstants.jl")

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

println("Arranging for eta_1 computation")
files_to_act_on = filter(contains(thickness_target), readdir("results"))

### Act to construct the epsilon_near_zero_dispersion
# Precompute the epsilon_near_zero_dispersion
wavevectors = LinRange(minimum_wavevector, maximum_wavevector, number_of_interpolation_bins);
initial_guess = [initial];
enz_frequencies = [ENZ.epsilon_near_zero_dispersion!(k, d, initial_guess) for (idx, k) in enumerate(wavevectors)]

ω_ENZ_int = linear_interpolation(wavevectors, enz_frequencies)

pol_frequencies = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[1] for k in wavevectors]
pol_frequencies = mapreduce(permutedims, vcat, pol_frequencies);

ω_p1_int = linear_interpolation(wavevectors, pol_frequencies[:, 1])
ω_p2_int = linear_interpolation(wavevectors, pol_frequencies[:, 2])
ω_p3_int = linear_interpolation(wavevectors, pol_frequencies[:, 3])

vg1(x) = Interpolations.gradient(ω_p1_int, x) / c_0
vg2(x) = Interpolations.gradient(ω_p2_int, x) / c_0
vg3(x) = Interpolations.gradient(ω_p3_int, x) / c_0

vgs(x) = [vg1(x), vg2(x), vg3(x)]

wavevectors = LinRange(minimum_wavevector, maximum_wavevector, number_of_plotting_bins);
wavenumbers = LinRange(minimum_wavenumber, maximum_wavenumber, number_of_plotting_bins);
frequencies = [wavenumber_to_omega(wn) for wn in wavenumbers]

pol_frequencies = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[1] for k in wavevectors]
@cast pol_frequencies[i,j] := pol_frequencies[i][j]
pol_eigenvectors = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[2] for k in wavevectors]
@cast pol_eigenvectors[i,j,k] := pol_eigenvectors[i][j,k]

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
    println(f"Evaluating eta_1 for thickness {d_term}nm at {t_term}K");

	save(f"figures/efficiency_{d_term}_{t_term}.pdf", current_figure())
end;
