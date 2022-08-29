import PhysicalConstants.CODATA2018: c_0

include("Constants.jl");

function omega_to_wavenumber(ω)
	k = ω / c_0;
	return k / (2. * π)
end

