module ENZ;

include("Materials.jl");
include("Constants.jl")
include("Conversions.jl");
include("TransferMatrix.jl");

using Optim;

function bilayer_sphp_dispersion(k_x)
	sqrt(
        c_0^2 * k_x^2 * (1. / ε_c + 1. / ε_∞) + ω_L^2
        - sqrt(
            (c_0^2 * k_x^2 * (1. / ε_c + 1. / ε_∞) + ω_L^2)^2
            - 4. * c_0^2 * k_x^2 * (ε_c * ω_T^2 + ε_∞ * ω_L^2) / ε_c
            / ε_∞
        )
    ) / sqrt(2)
end

function epsilon_near_zero_dispersion!(k, d, initial)
	thicknesses = [ustrip(d)];
	εs = [ε_Si, ε_SiC, ε_Si];

	F = optimize(ω->dispersion(ustrip(k), ω, thicknesses, εs), [initial[1]], Newton());

	result = Optim.minimizer(F)[1]
	initial[1] = result;
	return result * 1e14u"1/s"
end

function ω_n(k, d, m)
    n = 2 * m - 1
    return sqrt(
        ω_L^2 - β_L^2 * (k^2 + (n * π / d)^2)
        + 0.0u"1/s^2" * im
    )
end

function Ω_n(k, d, m, ω_ENZ)
    return 2 * β_L / d * sqrt(
        (ω_L^2 - ω_ENZ^2)
        / (ω_ENZ * ω_n(k, d, m))
    )
end

function bose_einstein(ω, T_lattice)
	1. / (exp(ħ * ω / (k_B * T_lattice)) - 1.)
end

function polariton_matrix!(k, d, initial, n_max)

	ω_ENZ = epsilon_near_zero_dispersion!(k, d, initial)

	polariton_matrix_from_enz_ω(k, d, ω_ENZ, n_max)
end


function polariton_matrix_from_enz_ω(k, d, ω_ENZ, n_max)
	matrix = zeros((n_max + 1, n_max + 1)) * 1u"1/s" * im

	ω_N = [ω_n(k, d, n) for n in 1:n_max];

	matrix[1, 1] = ω_ENZ;
	for (idx, ω) in enumerate(ω_N)
		matrix[idx+1, idx+1] = ω;
	end

	Ω_N = [Ω_n(k, d, n, ω_ENZ) for n in 1:n_max];
	for (idx, Ω) in enumerate(Ω_N)
		matrix[1, idx + 1] = im * Ω;
		matrix[idx + 1, 1] = - im * Ω;
	end

	matrix
end

function polariton_eigenvalues!(k, d, initial, n_max)
	matrix = polariton_matrix!(k, d, initial, n_max)
	eigs = eigen(ustrip(matrix));
	values = eigs.values;
	vectors = eigs.vectors
	return values * 1u"1/s", vectors
end

function polariton_eigenvalues_from_enz_ω(k, d, ω_ENZ, n_max)
	matrix = polariton_matrix_from_enz_ω(k, d, ω_ENZ, n_max)
	eigs = eigen(ustrip(matrix));
	values = eigs.values;
	vectors = eigs.vectors
	return values * 1u"1/s", vectors
end


function polariton_density_of_states(ω, k, d, initial, n_max)
	polariton_frequencies = polariton_eigenvalues!(k, d, initial, n_max)[1]
	dos = 0u"s"
	for i in 1:n_max+1
		dos += lorentzian_density_of_states(ω, polariton_frequencies[i], γ)
	end
	dos
end

function lorentzian_density_of_states(ω, ω_0, γ)
	return π * γ / 2. / ((ω - ω_0)^2 + (γ / 2)^2)
end

end