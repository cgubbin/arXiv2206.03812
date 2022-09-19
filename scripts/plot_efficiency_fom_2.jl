include("EvaluationConstants.jl")
include("../src/ENZ.jl")
include("../src/FieldBalance.jl")

using CairoMakie;
using DelimitedFiles
using .ENZ;
using .FieldBalance;
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

println("Arranging for computation of eta_2")
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
dk = wavevectors[2]-wavevectors[1];
wavenumbers = LinRange(minimum_wavenumber, maximum_wavenumber, number_of_plotting_bins);
frequencies = [wavenumber_to_omega(wn) for wn in wavenumbers]

pol_frequencies = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[1] for k in wavevectors]
@cast pol_frequencies[i,j] := pol_frequencies[i][j]
pol_eigenvectors = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[2] for k in wavevectors]
@cast pol_eigenvectors[i,j,k] := pol_eigenvectors[i][j,k]

n_thermal = [ENZ.bose_einstein(ω, T_lattice) for ω in pol_frequencies]

i = 1;
res = zeros((2, length(files_to_act_on)))
for (idx, file) in enumerate(files_to_act_on)
	mat = readdlm("results/" * file, ',', Float64);

	gamma_sum = sum(mat, dims=2)[1]



	substrings = split(split(file, ".")[1], "_")
	t_term = substrings[3]
	d_term = substrings[2]

    println(f"Evaluating eta_2 for thickness {d_term}nm at {t_term}K");

	T_e = parse(Float64, split(t_term, "K")[1]) * 1u"K"

	if T_e > T_lattice

		drift_velocity = FieldBalance.drift_velocity.(T_e, T_lattice);


		gamma_integrated = 2. * π * gamma_sum .* wavevectors / (2. * π)^2;
		result = sum(gamma_integrated) * dk * ENZ.γ / drift_velocity / 1e24u"1/m^3"

		res[1, idx] = ustrip(T_e); 
		res[2, idx] = upreferred(result)
	end


end


println("Generating plot of eta_2")
scatter(res[1, :], res[2, :]; color = "#56B4E9", linewidth = 2, label = L"\mathrm{Branch}\;1",
        axis = (xlabel = L"\mathrm{Electronic\;Temperature\;}(K)", ylabel = L"\mathrm{Figure\;of\;Merit\;\;}\eta^{(2)}", xgridcolor = :red,
            xlabelsize = 22, ylabelsize = 22,
            xgridstyle = :dashdot, xgridwidth = 0.85,
            xtickalign = 1, xticksize = 20),
        figure = (resolution = (600, 400), font = "CMU Serif"))

xlims!(300, 1900)
ylims!(0, 0.25)
save(f"figures/eta_2_{thickness_target}.pdf", current_figure())  
