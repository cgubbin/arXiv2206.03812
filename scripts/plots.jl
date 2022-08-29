include("../src/FieldBalance.jl")

using .FieldBalance;
using CairoMakie;
using Unitful;

lattice_temperature = 200u"K";

electronic_temperatures = LinRange(200, 2000, 50)u"K";

electric_fields = FieldBalance.electric_field.(electronic_temperatures, lattice_temperature);
drift_velocities = FieldBalance.drift_velocity.(electronic_temperatures, lattice_temperature);

fig = Figure(resolution=(1500, 500))
ax1, l1 = lines(fig[1, 1], ustrip(electronic_temperatures), ustrip(electric_fields), color = :red)
ax2, l2 = lines(fig[1, 2], ustrip(electronic_temperatures), ustrip(drift_velocities), color = :red)
ax3, l3 = lines(fig[1, 3], ustrip(electric_fields), ustrip(drift_velocities), color = :red)

linkxaxes!(ax1, ax2);
linkyaxes!(ax2, ax3);

xlims!(ax1, low=0, high=2000)
xlims!(ax2, low=0, high=2000)
xlims!(ax3, low=0)
ylims!(ax1, low=0)
ylims!(ax2, low=0)
ylims!(ax3, low=0)

ax1.xlabel = "Temperature (K)"
ax2.xlabel = "Temperature (K)"
ax3.xlabel = "Electric Field (V / m)"
ax1.ylabel = "Electric Field (V / m)"
ax2.ylabel = "Drift Velocity (m / s)"
ax3.ylabel = "Drift Velocity (m / s)"

fig

# d = 1e-9u"m";
# bulk_percentages = FieldBalance.bulk_like_absorption_percentage.(electronic_temperatures, lattice_temperature, d);
# layer_percentages = FieldBalance.layer_like_absorption_percentage.(electronic_temperatures, lattice_temperature, d);
# fig = Figure(resolution=(1500, 500))
# ax1, l1 = lines(fig[1, 1], ustrip(electronic_temperatures), bulk_percentages, color = :red)
# ax2, l2 = lines(fig[2, 1], ustrip(electronic_temperatures), layer_percentages, color = :red)
# fig