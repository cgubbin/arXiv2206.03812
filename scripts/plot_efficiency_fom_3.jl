include("../src/ENZ.jl")
include("../src/FieldBalance.jl")

using CairoMakie;
using DelimitedFiles
using .ENZ;
using .FieldBalance;
using Interpolations;
using LinearAlgebra;
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
dk = wavevectors[2]-wavevectors[1];
wavenumbers = LinRange(82500u"1/m", 100000u"1/m", n_bin_ω);
frequencies = [wavenumber_to_omega(wn) for wn in wavenumbers]

initial = [1.75];
pol_frequencies = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[1] for k in wavevectors]
@cast pol_frequencies[i,j] := pol_frequencies[i][j]
pol_eigenvectors = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[2] for k in wavevectors]
@cast pol_eigenvectors[i,j,k] := pol_eigenvectors[i][j,k]

vgs_d = [vgs(k) for k in wavevectors]
@cast vgs_d[i,j] := vgs_d[i][j]

T_lattice = 300u"K"

# n_thermal = [ENZ.bose_einstein(ω, T_lattice) for ω in pol_frequencies]

i = 1;
res = zeros((2, length(files_to_act_on)))
for (idx, file) in enumerate(files_to_act_on)
	mat = readdlm("results/" * file, ',', Float64);
	substrings = split(split(file, ".")[1], "_")
	t_term = substrings[3]
	d_term = substrings[2]

	T_e = parse(Float64, split(t_term, "K")[1]) * 1u"K"
	drift_velocity = FieldBalance.drift_velocity.(T_e, T_lattice);

	# gamma = ENZ.ħ * ENZ.γ .* vgs_d .* mat .* pol_frequencies;
	gamma = (ENZ.γ * mat) .* (ENZ.ħ * pol_frequencies) .* norm.(vgs_d);
	
	integral = 2. * π * gamma .* wavevectors / (2. * π)^2;
	result = sum(integral) * dk / drift_velocity / 1e24u"1/m^3"

	thermal_energy = 1.5 * ENZ.k_B * T_e;

	res[1, idx] = ustrip(T_e); 
	res[2, idx] = upreferred(result / thermal_energy)


end




scatter(res[1, :], res[2, :]; color = "#56B4E9", linewidth = 2, label = L"\mathrm{Branch}\;1",
        axis = (xlabel = L"\mathrm{Electronic\;Temperature\;}(K)", ylabel = L"\mathrm{Figure\;of\;Merit\;\;}\eta^{(3)}", xgridcolor = :red,
            xlabelsize = 22, ylabelsize = 22,
            xgridstyle = :dashdot, xgridwidth = 0.85,
            xtickalign = 1, xticksize = 20),
        figure = (resolution = (600, 400), font = "CMU Serif"))

xlims!(300, 1900)
ylims!(0, 1e-5)
save(f"figures/eta_3_{thickness_target}.pdf", current_figure())  