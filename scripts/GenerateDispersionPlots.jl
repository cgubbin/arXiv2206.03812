include("EvaluationConstants.jl")

using .ENZ;
using CairoMakie;
using FStrings;
using Interpolations;
using TensorCast;
using Unitful;

struct DispersionResult
	wavevectors::Vector{Float64}
	enz_frequencies::Vector{Float64}
	polariton_frequencies::Array{Float64, 2}
	polariton_eigenvectors::Array{ComplexF64, 3}
end

function generate_plot()
	# Generate a linear range of wavevectors over which to calculate the initial local dispersion
    println("Generating the ENZ interpolation")
	initial_guess = [initial];
	wavevectors = LinRange(minimum_wavevector, maximum_wavevector, number_of_interpolation_bins);
	enz_frequencies = [ENZ.epsilon_near_zero_dispersion!(k, d, initial_guess) for k in wavevectors]

	# Interpolate the dispersion to all intermediate wavevectors
	enz_interp = linear_interpolation(wavevectors, enz_frequencies)
	interp = [enz_interp(k) for k in wavevectors];

    println("Calculating the polariton frequencies")
	# Calculate the polariton frequencies
	initial_guess = [initial];
	pol_frequencies = [ENZ.polariton_eigenvalues!(k, d, initial_guess, n_max)[1] for k in wavevectors]
	@cast pol_frequencies[i, j] := pol_frequencies[i][j]
	pol_eigenvectors = [ENZ.polariton_eigenvalues!(k, d, initial_guess, n_max)[2] for k in wavevectors]
	@cast pol_eigenvectors[i, j, k] := pol_eigenvectors[i][j, k];

    println("Plotting the figure");
	# Plot limits and discretisation
	frequencies = LinRange(minimum_frequency, maximum_frequency, number_of_plotting_bins)
	wavevectors_map = LinRange(minimum_wavevector, maximum_wavevector, number_of_plotting_bins);
	initial_guess = [initial]
	density_of_states = [ENZ.polariton_density_of_states(ω, k, d, initial_guess, n_max) for k in wavevectors, ω in frequencies]
	density_of_states /= maximum(density_of_states); # Normalise -> for this plot we only care about the dispersion, not the values

	# t = "text";
	fig, _, hm = heatmap(
		ustrip(wavevectors_map) / 1e9, ustrip(ENZ.omega_to_wavenumber(frequencies)) / 100, ustrip(density_of_states),
		colormap = :lajolla, colorrange = (0, 1), interpolate = true,
		axis = (title = f"Dispersion for {thickness_nm}nm", xlabel = L"\mathrm{LTP\;Wavevector}\;(\mathrm{nm}^{-1})", ylabel = L"\mathrm{Wavenumber}\;(\mathrm{cm}^{-1})"),
		# axis = (title = L"Some %$(t) and some math: $\frac{2\alpha+1}{y}$"),
		font = "CMU Serif"
	);

	Colorbar(fig[:, end+1], hm)
	lines!(ustrip(wavevectors) / 1e9, ustrip(ENZ.omega_to_wavenumber(interp)) / 100, color = :green, linewidth=4)
	for i in 1:n_max+1
		lines!(ustrip(wavevectors) / 1e9, ustrip(ENZ.omega_to_wavenumber(pol_frequencies[:, i]))/ 100, color = :red, linewidth=4)
	end
	ylims!(low=825, high=1000)
	xlims!(low=ustrip(minimum_wavevector) / 1e9, high=ustrip(maximum_wavevector) / 1e9)

	fig, DispersionResult(ustrip(wavevectors), ustrip(enz_frequencies), ustrip(pol_frequencies), ustrip(pol_eigenvectors))
end

fig, results = generate_plot()
println("Saving dispersion plot to file")
save(f"figures/dispersion_{thickness_nm}nm.pdf", fig)
