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

	discrete_phonon_hopfield = discrete_polariton_hopfield[:, 2:end, :];

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

	interpolated_phonon_hopfields = [
		linear_interpolation(wavevectors, discrete_branch_phonon_hopfield)
		for discrete_branch_phonon_hopfields in eachslice(discrete_phonon_hopfield,dims=3) 
			for discrete_branch_phonon_hopfield in eachcol(discrete_branch_phonon_hopfields)
	]
	# interpolated_phonon_hopfields = transpose(reshape(
	# 	interpolated_phonon_hopfields, (n_max, n_max + 1)
	# )) # A matrix whose rows are the phonon hopfield coeffs for each branch
	# with shape (n_max + 1, n_max)


	interpolated_polariton_ω, interpolated_group_velocities, interpolated_enz_hopfield, interpolated_phonon_hopfields
end

polariton_ω, group_velocities, enz_hopfield, phonon_hopfields = generate_polariton_properties();
wavevectors = LinRange(minimum_wavevector, maximum_wavevector, number_of_plotting_bins);

# Branch 1
# println(size(phonon_hopfields), raw"size")
# println(phonon_hopfields, "phonon")
# println(size(enz_hopfield), "size")
# println(enz_hopfield, "enz")
# β1 = [abs(phonon_hopfields[1](q))^2 for q in wavevectors]
# β2 = [abs(phonon_hopfields[2](q))^2 for q in wavevectors]
# βe = [abs(enz_hopfield[1](q))^2 for q in wavevectors]
# println(β1 + β2 + βe, "1")

# q = wavevectors[20]
# j=1
# test = [abs(β(q))^2 for β in phonon_hopfields[(j-1)*n_max+1:j*n_max]]
# println(test, abs(enz_hopfield[j](q))^2)

# using PGFPlotsX

# figure = Plot(
# 	Table(ustrip(wavevectors) / 1e9, β1),
# 	Table(ustrip(wavevectors) / 1e9, β2),
# )
# pgfsave("./figures/test.pdf", figure)




if nprocs() == 1
	addprocs(40, topology=:master_worker, exeflags="--project=$(Base.active_project())")
end

@everywhere begin
  # instantiate environment
  println("Instantiating on each node");
  using Pkg;
  Pkg.instantiate();
  include("EvaluationConstants.jl")
  include("../src/RatesTyped.jl")

  using Interpolations;
  using Unitful;

  function generate_single_temperature(wavevector, temperature, polariton_ω, phonon_hopfields, group_velocities)
      result = [
                RatesTyped.gammaIntegratedkzj(wavevector, polariton_ω, phonon_hopfields, d, temperature, n_max, j)[1]
                for j in 1:n_max+1
               ];

	  result
  end

  function generate_single_temperature_absorption(wavevector, temperature, polariton_ω, phonon_hopfields, group_velocities)
      result = [
                RatesTyped.gammaIntegratedAbsorptionkzj(wavevector, polariton_ω, phonon_hopfields, d, temperature, n_max, j)[1]
                for j in 1:n_max+1
               ];

	  result
	end
end


println(f"Running on {nprocs()}")

for temperature in temperatures
	stripped_temperature = ustrip(temperature);
    println(f"Evaluating for temperature {stripped_temperature}K")
	result = @showprogress 1 f"Computing temperature {stripped_temperature}K..." pmap(wavevectors) do q
        try
            x = generate_single_temperature(
                q, temperature, polariton_ω, phonon_hopfields, group_velocities
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

    outfile = f"results/ep_{stripped_thickness:.0f}nm_{stripped_temperature:.0f}K.csv"
    writedlm(outfile, result, ",");

	result = @showprogress 1 f"Computing absorption for temperature {stripped_temperature}K..." pmap(wavevectors) do q
		try
			x = generate_single_temperature_absorption(
				q, temperature, polariton_ω, phonon_hopfields, group_velocities
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

	outfile = f"results/ep_abs_{stripped_thickness:.0f}nm_{stripped_temperature:.0f}K.csv"
	writedlm(outfile, result, ",");

end
