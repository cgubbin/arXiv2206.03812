include("EvaluationConstants.jl")
include("../src/RatesTyped.jl")

using CairoMakie;
using DelimitedFiles;
using .ENZ;
using Interpolations;
using TensorCast;
using Unitful;
using FStrings;
using ProgressMeter;
using ColorSchemes;

using Distributed

"""
	Generates the interpolated polariton frequencies and group group_velocities
"""
function generate_polariton_properties()
	wavevectors = LinRange(minimum_wavevector, maximum_wavevector, number_of_interpolation_bins);

	# Generate the interpolated epsilon near zero dispersion
	initial_guess = [initial];
	discrete_enz_ω = [ENZ.epsilon_near_zero_dispersion!(k, d, initial_guess) for (idx, k) in enumerate(wavevectors)];
	interpolated_enz_ω = linear_interpolation(wavevectors, discrete_enz_ω);

	# Generate the discrete polariton frequencies
	discrete_polariton_ω = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, interpolated_enz_ω(k), n_max)[1] for k in wavevectors]
	@cast discrete_polariton_ω[i, j] := discrete_polariton_ω[i][j];

	discrete_polariton_hopfield = [ENZ.polariton_eigenvalues_from_enz_ω(k, d, interpolated_enz_ω(k), n_max)[2] for k in wavevectors]
	@cast discrete_polariton_hopfield[i, j, k] := discrete_polariton_hopfield[i][j, k];
	discrete_enz_hopfield = discrete_polariton_hopfield[:, 1, :];

	# Form an array of interpolated polariton frequencies
	interpolated_polariton_ω = [
		linear_interpolation(wavevectors, discrete_branch_ω) 
		for discrete_branch_ω in eachcol(discrete_polariton_ω)
			];

	interpolated_group_velocities = [
		x->Interpolations.gradient(interpolated_branch_ω, x)
		for interpolated_branch_ω in interpolated_polariton_ω
	]	;

	interpolated_enz_hopfield = [
		linear_interpolation(wavevectors, discrete_branch_enz_hopfield) 
		for discrete_branch_enz_hopfield in eachcol(discrete_enz_hopfield)
	];


	interpolated_polariton_ω, interpolated_group_velocities, interpolated_enz_hopfield
end

polariton_ω, group_velocities, enz_hopfield = generate_polariton_properties();
wavevectors = LinRange(minimum_wavevector, maximum_wavevector, number_of_plotting_bins);

if nprocs() == 1
	addprocs(8, topology=:master_worker, exeflags="--project=$(Base.active_project())")
end

@everywhere begin
  # instantiate environment
  using Pkg; Pkg.instantiate()
	include("EvaluationConstants.jl")
	include("../src/RatesTyped.jl")

	using Interpolations;
	using Unitful;

	function generate_single_temperature(wavevector, temperature, polariton_ω, enz_hopfield, group_velocities)
		result = [
			RatesTyped.gammaIntegratedkzj(wavevector, polariton_ω, enz_hopfield, d, temperature, n_max, j)[1]
			for j in 1:n_max+1
		];

		result
	end

end


println(nprocs())
for temperature in temperatures
	stripped_temperature = ustrip(temperature);
	result = @showprogress 1 f"Computing temperature {stripped_temperature}K..." pmap(wavevectors) do q
		try
			x = generate_single_temperature(
				q, temperature, polariton_ω, enz_hopfield, group_velocities
			);
			x
		catch e
			@warn "failed to solve"
			rethrow(e)
		end
	end

	@cast result[i, j] := result[i][j];

	result = result * 1e24u"1/s" / ENZ.γ;

	stripped_thickness = ustrip(d) / 1e-9;

	outfile = f"results_tmp/ep_{stripped_thickness:.0f}nm_{stripped_temperature:.0f}K.csv"
	writedlm(outfile, result, ",");

end
