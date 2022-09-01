include("../src/ENZ.jl")

using Unitful;
import PhysicalConstants.CODATA2018: c_0;

function wavenumber_to_omega(wn)
	k = wn * 2. * π;
	return k * c_0
end

const thickness_nm = 3; # The target thickness for the simulation in nm
const d = thickness_nm * 1e-9u"m" # The target thickness for the simulation
const n_max = 3; # The number of phonon modes to include in the problem
const minimum_wavenumber = 82500u"1/m" # The minumum wavenumber to use in plotting
const maximum_wavenumber = 1e5u"1/m" # The maximum wavenumber to use in plotting
const minimum_wavevector = ENZ.ω_L * ENZ.ε_c / c_0 * 1.4; # The minimum wavevector (skipping the invalid low-wavevector rise)
const maximum_wavevector = 1e9u"1/m"; # The cut-off wavevector for plotting and calculations
const minimum_frequency = wavenumber_to_omega(minimum_wavenumber);
const maximum_frequency = wavenumber_to_omega(maximum_wavenumber);
const number_of_interpolation_bins = 1000; # The number of bins for constructing the interpolation functions
const number_of_plotting_bins = 200; # The number of bins for plot output
const initial = 1.75 # The inital frequency guess in 100THz, for calculating the mode dispersion
const T_lattice = 300u"K"; # The lattice temperature in Kelvin

const temperatures = LinRange(400, 2000, 17) * 1u"K"; # The default temperatures to run the calculation over
