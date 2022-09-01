import PhysicalConstants.CODATA2018: c_0

include("Constants.jl");

"""
	omega_to_wavenumber(ω)

	Converts the ω, a frequency in "rad/s" to a wavenumber in inverse centimeters
"""
function omega_to_wavenumber(ω)
	k = ω / c_0;
	return k / (2. * π)
end

