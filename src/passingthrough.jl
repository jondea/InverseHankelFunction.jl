
struct PassingThrough{T}
    point::T
end

import Base.isfinite
isfinite(p::PassingThrough) = isfinite(p.point)

import Base.iszero
iszero(p::PassingThrough) = iszero(p.point)

import Base.isnan
isnan(p::PassingThrough) = isnan(p.point)

"""
    invhankelh1(ν, h, p::PassingThrough)

Find z such that H^{(1)}_\\nu(z) = h and the solution branch passes through p, the inverse of SpecialFunctions.hankelh1
"""
function invhankelh1(ν::Integer, h::Number, p::PassingThrough)

    # Handle some special cases
    if isnan(p)
        throw(DomainError(p.point, "invhankelh1 never passes through $(p.point) because it is not a number"))
    end
    if !isfinite(p)
        throw(DomainError(p.point, "Multiple branches of invhankelh1 pass through non finite ($(p.point)) find another value that this branch passes through"))
    end
    if iszero(p)
        throw(DomainError(p.point, "Multiple branches of invhankelh1 pass through zero ($(p.point)) find another value that this branch passes through"))
    end

    # Starting point for our continuation
    z₀ = p.point

    return invhankelh1n(ν, z₀, h/hankelh1(ν,z₀))
end
