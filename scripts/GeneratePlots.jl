include("../src/ENZ.jl")
include("EvaluationConstants.jl")
include("../src/FieldBalance.jl")

using CairoMakie;
using DelimitedFiles
using .ENZ;
using Interpolations;
using .FieldBalance;
using TensorCast;
using Unitful;
using FStrings;
using ColorSchemes;

import PhysicalConstants.CODATA2018: c_0;

function wavenumber_to_omega(wn)
	k = wn * 2. * π;
	return k * c_0
end

function arrange()
	println("Arranging for plots computation")
	files_to_act_on = filter(contains(thickness_target), readdir("results"))

	### Act to construct the epsilon_near_zero_dispersion
	# Precompute the epsilon_near_zero_dispersion
	wavevectors = LinRange(minimum_wavevector, maximum_wavevector, number_of_interpolation_bins);
	initial_guess = [initial];
	enz_frequencies = [ENZ.epsilon_near_zero_dispersion!(k, d, initial_guess) for (idx, k) in enumerate(wavevectors)]

	ω_ENZ_int = linear_interpolation(wavevectors, enz_frequencies)

	pol_frequencies = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[1] for k in wavevectors]
	pol_frequencies = mapreduce(permutedims, vcat, pol_frequencies);

	ω_int = [linear_interpolation(wavevectors, pol_frequencies[:, i]) for i in 1:n_max + 1]
	vgs(x) = [Interpolations.gradient(ω_int[i], x) / c_0 for i in 1:n_max + 1]

	wavevectors = LinRange(minimum_wavevector, maximum_wavevector, number_of_plotting_bins);
	wavenumbers = LinRange(minimum_wavenumber, maximum_wavenumber, number_of_plotting_bins);

	vgs_d = [vgs(k) for k in wavevectors]
	@cast vgs_d[i,j] := vgs_d[i][j]

	pol_frequencies = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[1] for k in wavevectors]
	@cast pol_frequencies[i,j] := pol_frequencies[i][j]
	pol_eigenvectors = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ_int(k), n_max)[2] for k in wavevectors]
	@cast pol_eigenvectors[i,j,k] := pol_eigenvectors[i][j,k]

	n_thermal = [ENZ.bose_einstein(ω, T_lattice) for ω in pol_frequencies]

	return wavevectors, wavenumbers, pol_frequencies, vgs_d, n_thermal
end

function generate_eta_1_plot(wavevectors, n_thermal)
	files_to_act_on = filter(contains(thickness_target), readdir("results"))
	for file in files_to_act_on
		mat = readdlm("results/" * file, ',', Float64);
		lines(ustrip(wavevectors), mat[:, 1] ./ n_thermal[:, 1]; color = "#56B4E9", linewidth = 2, label = L"\mathrm{Branch}\;1",
					axis = (xlabel = L"\mathrm{LTP\;Wavevector}\;(\mathrm{m}^{-1})", ylabel = L"\mathrm{Figure\;of\;Merit\;\;}\eta^{(1)}", xgridcolor = :red,
							xlabelsize = 22, ylabelsize = 22,
							xgridstyle = :dashdot, xgridwidth = 0.85,
							xtickalign = 1, xticksize = 20),
					figure = (resolution = (600, 400), font = "CMU Serif"))
		# break
		for i in 2:n_max
			lines!(ustrip(wavevectors), mat[:, i] ./ n_thermal[:, i]; color = :black, linestyle = :dash, label = L"\mathrm{Branch}\;2")
		end
		xlims!(minimum(ustrip(wavevectors)), maximum(ustrip(wavevectors)))
		# ylims!(0, 50)
		axislegend("Legend", position = :lt)
		substrings = split(split(file, ".")[1], "_")
		t_term = substrings[3]
		d_term = substrings[2]
		println(f"Evaluating eta_1 for thickness {d_term} at {t_term}");

		outfile = f"postprocessing/eta1_{d_term}_{t_term}.csv"
		writedlm(outfile, mat ./ n_thermal, ",");
	
		save(f"figures/efficiency_{d_term}_{t_term}.pdf", current_figure())
	end;
	
	outfile = f"postprocessing/eta1_wavevector.csv"
	writedlm(outfile, ustrip(wavevectors) / 1e9, ",");
end

function generate_eta_2_plot(wavevectors)
	files_to_act_on = filter(contains(thickness_target), readdir("results"))
	res = zeros((2, length(files_to_act_on)))

	dk = wavevectors[2] - wavevectors[1];

	for (idx, file) in enumerate(files_to_act_on)
		mat = readdlm("results/" * file, ',', Float64);

		gamma_sum = sum(mat, dims=2)[1]



		substrings = split(split(file, ".")[1], "_")
		t_term = substrings[3]
		d_term = substrings[2]

		println(f"Evaluating eta_2 for thickness {d_term} at {t_term}");

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

	# xlims!(300, 1900)
	# ylims!(0, 0.25)
	save(f"figures/eta_2_{thickness_target}.pdf", current_figure())  

	outfile = f"postprocessing/eta2.csv"
	writedlm(outfile, res, ",");
end

function generate_eta_3_plot(wavevectors, pol_frequencies, vgs_d)
	files_to_act_on = filter(contains(thickness_target), readdir("results"))
	res = zeros((2, length(files_to_act_on)))
	dk = wavevectors[2] - wavevectors[1];

	for (idx, file) in enumerate(files_to_act_on)
		mat = readdlm("results/" * file, ',', Float64);
		substrings = split(split(file, ".")[1], "_")
		t_term = substrings[3]
		d_term = substrings[2]
	
		println(f"Evaluating eta_3 for thickness {d_term} at {t_term}");
	
		T_e = parse(Float64, split(t_term, "K")[1]) * 1u"K"
		drift_velocity = FieldBalance.drift_velocity.(T_e, T_lattice);
	
		gamma = - (ENZ.γ * mat) .* (ENZ.ħ * pol_frequencies) .* vgs_d;
		
		integral = 2. * π * gamma .* wavevectors / (2. * π)^2;
		result = sum(integral) * dk / drift_velocity / 1e24u"1/m^3"
	
		thermal_energy = 1.5 * ENZ.k_B * T_e;
	
		res[1, idx] = ustrip(T_e); 
		res[2, idx] = upreferred(result[1] / thermal_energy)
	
	
	end

	println("Generating plot of eta_3")
	
	scatter(res[1, :], res[2, :]; color = "#56B4E9", linewidth = 2, label = L"\mathrm{Branch}\;1",
					axis = (xlabel = L"\mathrm{Electronic\;Temperature\;}(K)", ylabel = L"\mathrm{Figure\;of\;Merit\;\;}\eta^{(3)}", xgridcolor = :red,
							xlabelsize = 22, ylabelsize = 22,
							xgridstyle = :dashdot, xgridwidth = 0.85,
							xtickalign = 1, xticksize = 20),
					figure = (resolution = (600, 400), font = "CMU Serif"))
	
	xlims!(300, 1900)
	# ylims!(0, 1e-5)
	save(f"figures/eta_3_{thickness_target}.pdf", current_figure())  

	outfile = f"postprocessing/eta3.csv"
	writedlm(outfile, res, ",");
	
end

function generate_figure_4(wavevectors, wavenumbers, pol_frequencies, vgs_d)
	frequencies = [wavenumber_to_omega(wn) for wn in wavenumbers]
	files_to_act_on = filter(contains(thickness_target), readdir("results"))
	for file in files_to_act_on
		println(file)
		mat = readdlm("results/" * file, ',', Float64);
	
		electrical_map = [
			sum([mat[j, m] * ENZ.lorentzian_density_of_states(ω, pol_frequencies[j, m], ENZ.γ) for m in 1:n_max+1])
			for j in 1:number_of_plotting_bins, ω in frequencies 
		]
	
		electrical_map_vg = [
			sum([abs(vgs_d[j, m][1]) * mat[j, m] * ENZ.lorentzian_density_of_states(ω, pol_frequencies[j, m], ENZ.γ) for m in 1:n_max+1])
			for j in 1:number_of_plotting_bins, ω in frequencies 
		]
	
		thermal_map = [
			sum([ENZ.bose_einstein(pol_frequencies[j, m], T_lattice) * ENZ.lorentzian_density_of_states(ω, pol_frequencies[j, m], ENZ.γ) for m in 1:n_max+1])
			for j in 1:number_of_plotting_bins, ω in frequencies 
		]
	
		thermal_map_vg = [
			sum([abs(vgs_d[j, m][1]) * ENZ.bose_einstein(pol_frequencies[j, m], T_lattice) * ENZ.lorentzian_density_of_states(ω, pol_frequencies[j, m], ENZ.γ) for m in 1:n_max+1])
			for j in 1:number_of_plotting_bins, ω in frequencies 
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

		outfile = f"postprocessing/figure4_a_map_{d_term}_{t_term}.csv"
		writedlm(outfile, ustrip(thermal_map), ",");
		outfile = f"postprocessing/figure4_b_map_{d_term}_{t_term}.csv"
		writedlm(outfile, ustrip(thermal_map_vg), ",");
		outfile = f"postprocessing/figure4_c_map_{d_term}_{t_term}.csv"
		writedlm(outfile, ustrip(electrical_map), ",");
		outfile = f"postprocessing/figure4_d_map_{d_term}_{t_term}.csv"
		writedlm(outfile, ustrip(electrical_map_vg), ",");
	end;	
	outfile = f"postprocessing/figure4_x.csv"
	writedlm(outfile, ustrip(wavevectors) / 1e9, ",");
	outfile = f"postprocessing/figure4_y.csv"
	writedlm(outfile, ustrip(wavenumbers) / 100, ",");
end

println("Arranging")
wavevectors, wavenumbers, pol_frequencies, vgs_d, n_thermal = arrange();
println("Plotting eta 1")
generate_eta_1_plot(wavevectors, n_thermal)
println("Plotting eta 2")
generate_eta_2_plot(wavevectors)
println("Plotting eta 3")
generate_eta_3_plot(wavevectors, pol_frequencies, vgs_d)
println("Plotting figure 4")
generate_figure_4(wavevectors, wavenumbers, pol_frequencies, vgs_d)