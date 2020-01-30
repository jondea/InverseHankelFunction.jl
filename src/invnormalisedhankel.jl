
"""
    normalisedhankelh1(ν, z₀ [,z])

Hankel function normalised at some z₀, H^{(1)}_\\nu(z)/H^{(1)}_\\nu(z₀), returns a function if `z` is not provided
"""
function normalisedhankelh1(ν, z₀, z)
    hankelh1(ν, z) / hankeh1(ν, z₀)
end

"""
Struct to cache the value of the hankel function at `z₀`, can be called like a function after construction
"""
struct NormalisedHankelH1
    ν::Complex{Float64}
    h₀::Complex{Float64}
end

normalisedhankelh1(ν, z₀) = NormalisedHankelH1(ν, hankeh1(ν, z₀))

(nh::NormalisedHankelH1)(z) = hankelh1(nh.ν, z)/nh.h₀

# "The gradient of the field as an ODE"
# dtν_dt(field::Field) = - hankel(nh.ν, z) / diffhankelh1(nh.ν, z)
#
# dtν_dν(nh::NormalisedHankelH1, z)::Complex = - nh.h₀ / diffhankelh1(nh.ν, z)
#
# dtν_dξ(nh::NormalisedHankelH1, z)::Complex = nh.h₀ / diffhankelh1(nh.ν, z)

@doc raw"""
    invnormalisedhankelh1(ν::Integer, h̄::Number, z₀::Number)

Find `z` such that ``H^{(1)}_\\nu(z)/H^{(1)}_\\nu(z₀) = \bar{h}`` for the branch continued from `z₀`
See also: [`normalisedhankelh1`](@ref)
"""
function invnormalisedhankelh1(ν, hbar::Number, z₀::Number)

    # SMALL_ARGUMENT_THRESHOLD = 0.5
    # LARGE_ARGUMENT_THRESHOLD = 2
    # if hankel_arg_asymptotic_scale(ν, z₀) < SMALL_ARGUMENT_THRESHOLD
    #     # Work out this branch from the nearest root of unity of h
    #     z = inv_small_arg_hankelh1(ν, h, )
    #     # Only works as far asymptotic scale for z is still below threshold
    # elseif hankel_arg_asymptotic_scale(ν, z₀) > LARGE_ARGUMENT_THRESHOLD
    #     # Work out branch from nearest 2πb to z₀
    #     z = inv_large_arg_hankelh1(ν, h, )
    #     # Only works as far asymptotic scale for z is still above threshold
    # end

    # Now do some numerical continuation until we have reached z, or we are back
    # In a large/small argument asymptotic regime

    return invnormalisedhankel_adaptive_solve(ν, z₀, hbar; householder_order=2, ε=1.0e-15,
        show_trace=false, step_max=0.1, ζ_jump_max=0.5, dζ_dξ_angle_jump_max=π/8, N_iter_max=10, silent_failure=false)
end

# vector version which reuses
function invnormalisedhankelh1_sortedvec(ν, hbars::AbstractVector{<:Real}, z₀::Number)
    zs = similar(hbars, complex(eltype(hbars)))
    z = z₀
    hbar_prev = one(eltype(hbars))
    for (i,hbar) in enumerate(hbars)
        z = invnormalisedhankel_adaptive_solve(ν, z₀, hbar, z, hbar_prev)
        zs[i] = z
        hbar_prev = hbar
    end
    return zs
end

# function invnormalisedhankel_adaptive_solve(ν::Number, z₀::Number, ξ::Number)
#     nh_fnc = normalisedhankelh1(ν, z₀)
#     invnormalisedhankel_adaptive_solve(nh_fnc::Function, ξ::Number)
# end

function invnormalisedhankel_adaptive_solve(ν::Number, z₀::Number, ξ_target::Number, z::Number=z₀, ξ::Number=one(ξ_target); householder_order=2, ε=1.0e-12,
    show_trace=false, step_max=0.1, ζ_jump_max=0.5, dζ_dξ_angle_jump_max=π/8, N_iter_max=10, silent_failure=false)

    h₀ = hankelh1(ν, z₀)

    # Last resort kill switch
    overall_iter = 0

    # Start at z and ξ
    ζ = z - z₀

    h = hankelh1(ν, z)
    h_ν_minus_1 = hankelh1(ν-1, z)
    dh_dζ = h_ν_minus_1 - ν/z*h

    # Tangent of ζ/z
    dζ_dξ = h₀/dh_dζ

    dζ_dξ_prev = dζ_dξ

    ζ_prev = ζ
    ξ_prev = ξ

    # Try to guess what will be a good size for the next step without going past our max
    ξ_step = sign(ξ_target-ξ)*min(step_max, 0.9*ζ_jump_max/abs.(dh_dζ), abs(ξ_target-ξ))

    # Keep going until we get to our target
    while abs(ξ_target-ξ) > 0

        # Last resort kill switch
        overall_iter += 1
        if overall_iter > 1000
            if silent_failure
                break
            else
                error("Too many steps overall, quitting")
            end
        end

        # Step in ξ
        ξ = ξ_prev + ξ_step

        # Predictor step in ζ or z using tangent (Euler)
        # ζ = ζ_prev + ξ_step*dζ_dξ

        ζ = ζ_prev

        # Calculate new point (z) and hankel function/its derivatives at new point
        z = z₀ + ζ
        h = hankelh1(ν, z)
        h_ν_minus_1 = hankelh1(ν-1, z)
        dh_dζ = h_ν_minus_1 - ν/z*h

        residual = h/h₀ - ξ

        # Corrector iterations
        iter = 0
        while abs(residual) > ε
            dresidual_dζ = dh_dζ/h₀
            δ = - residual/dresidual_dζ


            if householder_order == 1
                ζ = ζ + δ
            elseif householder_order == 2
                h_ν_minus_2 = hankelh1(ν-2, z)
                dh_ν_minus_1_dζ = h_ν_minus_2 - (ν-1)/z*h_ν_minus_1
                d2residual_dζ2 = dh_ν_minus_1_dζ -ν/z*dh_dζ + ν/z^2*h
                ζ = ζ + δ*(1 + δ*(d2residual_dζ2/(2*dresidual_dζ)))^-1
            elseif householder_order == 3
                h_ν_minus_2 = hankelh1(ν-2, z)
                dh_ν_minus_1_dζ = h_ν_minus_2 - (ν-1)/z*h_ν_minus_1
                d2residual_dζ2 = dh_ν_minus_1_dζ -ν/z*dh_dζ + ν/z^2*h
                error("Implement d3residual_dζ3")
                ζ = ζ + δ*(1 + δ*(d2residual_dζ2/(2*dresidual_dζ))) / (1 + (d2residual_dζ2/dresidual_dζ)*δ + 1//6*(d3residual_dζ3/dresidual_dζ)*δ^2)
            else
                error("Householder order $householder_order not implemented. Use 1, 2 or 3.")
            end

            # Calculate new point (z) and hankel function/its derivatives at new point
            z = z₀ + ζ
            h = hankelh1(ν, z)
            h_ν_minus_1 = hankelh1(ν-1, z)
            dh_dζ = h_ν_minus_1 - ν/z*h

            residual = h/h₀ - ξ

            # Keep track of number of iterations and exit to reduce stepsize if we do too many
            iter += 1
            if iter > N_iter_max
                break
            end
        end

        # Get tangent of ζ/z at new point
        dζ_dξ = h₀/dh_dζ

        # Get angle between this and previous
        # Angle between a=re^(im θ) and b=ρe^(im φ) is angle(r/ρ e^(im (θ-φ)))
        dζ_dξ_angle_jump = abs(angle(dζ_dξ_prev/dζ_dξ))
        residual = h/h₀ - ξ

        # Only accept this stepsize if the tangent and value hasn't changed too much
        if (   dζ_dξ_angle_jump <= dζ_dξ_angle_jump_max
            && abs(ζ-ζ_prev) <= ζ_jump_max
            && abs(residual) <= ε )

            # Reset stepsize, try to guess what will be a good size for the next step without going past our max
            ξ_step = sign(ξ_target-ξ)*min(step_max, 0.9*ζ_jump_max/abs.(dh_dζ), abs(ξ_target-ξ))

            # Store current as the starting point for next step
            dζ_dξ_prev = dζ_dξ
            ζ_prev = ζ
            ξ_prev = ξ
        else

            # Reject the step and reduce step size
            # Go back to previous
            dζ_dξ = dζ_dξ_prev
            ζ = ζ_prev
            ξ = ξ_prev

            # Half the size of ξ_step
            ξ_step = ξ_step/2
        end

    end

    return z
end
