function ε_SiC(ω)
    return ε_Lorentz(ω, ω_L, ω_T, ε_∞, γ)
end

function ε_Si(wavenumber)
    return ε_c
end

function ε_Lorentz(ω, ω_L, ω_T, ε_∞, γ)
    return ε_∞ * (
        (ω_L^2 - ω * (ω - im * γ))
        / (ω_T^2 - ω * (ω - im * γ))
    )
end
