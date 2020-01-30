
"""
    hankel_arg_asymptotic_scale(ν, z)

Determines largeness or smallness of argument to determine the validity of the asymptotic forms of the hankel functions.
"""
function hankel_arg_asymptotic_scale(ν, z)
    return abs(z)/sqrt(abs(ν)+1)
end

"""
    small_arg_hankelh1(ν::Integer, z::Number) -> Complex

The small argument (`z`) asymptotic of `hankelh1`
"""
function small_arg_hankelh1(ν::Integer, z::Number)
# function f(ν::Integer, z::Number)
    if iszero(ν)
        # Range: whole imaginary line, real part >= 1-2 and =< 1+2
        return  1 + 2im/π * log(z)
    else
        # Range: whole complex plane
        # May need to compute second two terms together as a product for numerical stability
        return - im/π * gamma(ν) * (z/2)^(-abs(ν))
    end
end


"""
    inv_small_arg_hankelh1(ν::Integer, h::Complex, b::Integer) -> Complex

The `b`th branch of the inverse of the small argument (`z`) asymptotic of `hankelh1`
Not valid for b > ν
"""
function inv_small_arg_hankelh1(ν::Integer, h::Complex{T}, b::Integer) where T
    if iszero(ν)
        if (real(h) > 3 || real(h) < -1) # Outside of the range of the function we are inverting, therefore a domain error
            throw(DomainError(h, "For inv_small_arg_hankelh1(0, h, b), -1 < real(h) < 3"))
        end
        return exp((h-1)*π/2im)
    else
        # Branch chooses root of unity
        return  2 * (h/(-im/π*gamma(ν)))^(-1/abs(ν)) * exp(2π * im * -b/abs(ν))
    end
end


"""
    large_arg_hankelh1(ν::Integer, z::Number) -> Complex

The large argument (`z`) asymptotic of `hankelh1`
"""
function large_arg_hankelh1(ν::Integer, z::Number)
    sqrt(2/(π*z)) * exp(im*(z - (ν*π)/2 - π/4))
end

"""
    inv_large_arg_hankelh1(ν::Integer, h::Complex, b::Integer) -> Complex

The `b`th branch of the inverse of the large argument (`z`) asymptotic of `hankelh1`
"""
function inv_large_arg_hankelh1(ν::Integer, h::Number, b::Integer)
    # Prefactor in Hankel asymptotic
    C = sqrt(2/π) * exp(im*(- (ν*π)/2 - π/4))
    # Even branches of LambertW don't seem to work, use 2b+1, investigate further
    return (im/2) * lambertw(-2im * (h/C)^(-2), 2b+1)
end

# Version using logs, is the first order approximation of LambertW form
# "Inverse of hankel function using large argument asymptotic form"
# function inv_large_arg_hankelh1(ν::Integer, h::Number, b::Integer)
#     return -im*log(h) + 2π*b + (ν*π)/2 + π/4 # + log log terms? + o(1)
# end
