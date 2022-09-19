using Pkg;
Pkg.instantiate();

println("Running dispersion generator")
include("GenerateDispersionPlots.jl")

println("Generating the Gamma")
include("GenerateGamma.jl")

println("Plotting the figures")
include("GeneratePlots.jl")
