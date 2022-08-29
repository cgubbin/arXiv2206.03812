import PhysicalConstants.CODATA2018: c_0;
using LinearAlgebra;

function κ(k, ω, ε)
	sqrt(ε(ω) * ω^2 / c_0^2 - k^2 + 0.0u"1/m^2" * im)
end

function λ(k, ω, ε)
	sqrt(k^2 - ε(ω) * ω^2 / c_0^2 + 0.0u"1/m^2" * im)
end


function generateM(k, ω, ε, d)
	M_i = ε(ω);
	κ_i = κ(k, ω, ε);
	κ_id = κ_i * d;
	M_11 = cos(κ_id);
	M_12 = -sin(κ_id) * M_i / κ_i;
	M_21 = sin(κ_id) * κ_i / M_i;
	M_22 = cos(κ_id);
	return [M_11 M_12; M_21 M_22]
end

function dispersion(k, ω, vd, vε)
	k = k * 1u"1/m";

	ω = ω[1] * 1e14u"1/s";
	ε = vε[1];
	λ_1 = λ(k, ω, ε);
	ε_1 = ε(ω);
	ϕ = [1.0; -λ_1 / ε_1];

	for (ε, d) in zip(vε[2:end], vd)
		d = d * 1u"m";
		M = generateM(k, ω, ε, d);
		ϕ = M * ϕ;
	end

	ε = vε[end];
	λ_n = λ(k, ω, ε);
	ε_n = ε(ω);
	F = abs(ϕ[1] - ϕ[2] * ε_n / λ_n)
	return F
end
