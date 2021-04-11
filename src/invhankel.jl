#=

"""
    invhankelh1(ν, h, b::Integer)

Find z such that H^{(1)}_\\nu(z) = h for the branch indexed by b, this is the inverse of SpecialFunctions.hankelh1
"""
function invhankelh1(ν::Integer, h::Number, b::Integer)

    z = complex(zero(h))

    SMALL_ARGUMENT_THRESHOLD = 0.5
    LARGE_ARGUMENT_THRESHOLD = 2
    if hankel_arg_asymptotic_scale(ν, z) < SMALL_ARGUMENT_THRESHOLD
        z = inv_small_arg_hankelh1(ν, h, b)
    elseif hankel_arg_asymptotic_scale(ν, z) > LARGE_ARGUMENT_THRESHOLD
        z = inv_large_arg_hankelh1()
    end

end

# Should we override broadcast to make it serial and reuse close values
=#
