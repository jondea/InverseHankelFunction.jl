"""
    hankelh1n(ν, z₀ [,z])

Hankel function normalised at some z₀, H^{(1)}_\\nu(z)/H^{(1)}_\\nu(z₀), returns a function if `z` is not provided
"""
hankelh1n(ν, z₀, z) = hankelh1(ν, z) / hankelh1(ν, z₀)

"""
Struct to cache the value of the Hankel function of the first kind at `z₀`, can be called like a function after construction
"""
struct HankelH1N
    "Order of the normalised Hankel function"
    ν::Float64
    "Normalise the Hankel function everywhere by the value of the Hankel function at this point in complex space"
    z₀::Complex{Float64}
    "Hankel function evaluated at z₀, we use this value to normalise the Hankel function"
    h₀::Complex{Float64}
end

HankelH1N(ν, z₀) = HankelH1N(ν, z₀, hankelh1(ν, z₀))

hankelh1n(ν, z₀) = HankelH1N(ν, z₀)

(nh::HankelH1N)(z) = hankelh1(nh.ν, z)/nh.h₀
