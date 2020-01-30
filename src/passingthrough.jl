
struct PassingThrough{T}
    point::T
end

isfinite(p::PassingThrough) = isfinite(p.point)
iszero(p::PassingThrough) = iszero(p.point)


"""
    invhankelh1(ν, h, p::PassingThrough)

Find z such that H^{(1)}_\\nu(z) = h and the solution branch passes through p, the inverse of SpecialFunctions.hankelh1
"""
function invhankelh1(ν::Integer, h::Number, p::PassingThrough)

    # Handle some special cases
    if isinf(p)
        throw(DomainError(p.point, "Multiple branches of invhankelh1 pass through $(p.point), find another value that this branch passes through"))
    end
    if iszero(p)
        throw(DomainError(p.point, "Multiple branches of invhankelh1 pass through $(p.point), find another value that this branch passes through"))
    end
    if isnan(p)
        throw(DomainError(p.point, "invhankelh1 never passes through $(p.point)"))
    end

    # Starting point for our continuation
    z = p.point

    SMALL_ARGUMENT_THRESHOLD = 0.5
    LARGE_ARGUMENT_THRESHOLD = 2
    if abs(z)/sqrt(abs(ν)+1) < SMALL_ARGUMENT_THRESHOLD

    elseif abs(z)/sqrt(abs(ν)+1) > LARGE_ARGUMENT_THRESHOLD

    end



    return 1.0

end
