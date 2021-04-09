using InverseHankelFunction
import Test: @test, @testset, @test_throws
import SpecialFunctions: hankelh1
import BenchmarkTools: @benchmark

include("utils.jl")

@testset "Hankel function asymptotics       " begin

    @testset "Small argument" begin

        # Some small values to test our small argument asymptotics
        small_values = 1e-4 .*[1+0im, 1im, -1+0im, -1im]

        @testset "H ≈ H_asymptotic" begin
            @test small_arg_hankelh1.(0, small_values)  ≈ hankelh1.(0, small_values)  rtol=0.1
            @test small_arg_hankelh1.(1, small_values)  ≈ hankelh1.(1, small_values)  rtol=1e-5
            @test small_arg_hankelh1.(2, small_values)  ≈ hankelh1.(2, small_values)  rtol=1e-5
            @test small_arg_hankelh1.(4, small_values)  ≈ hankelh1.(4, small_values)  rtol=1e-5
            @test small_arg_hankelh1.(8, small_values)  ≈ hankelh1.(8, small_values)  rtol=1e-5
            @test small_arg_hankelh1.(11, small_values) ≈ hankelh1.(11, small_values) rtol=1e-5
        end

        # Test whether f f⁻¹ = identity for different Hankel order (n) and solution branches (b)
        # We apply f⁻¹ first because action of f elliminates effect of branch
        @testset "f f⁻¹(z) = z" begin
            for n in [0,1,2,4,8,11]
                for b in [0,1,2,3,4,5,10,11]
                    f   = z->   small_arg_hankelh1(n, z)
                    invf = z->inv_small_arg_hankelh1(n, z, b)
                    # f is small argument asymptotic, so f(small_values) should be in domain of invf
                    testvec = f.(small_values)
                    if !all(isapprox.((f  ∘ invf).(testvec), testvec; rtol=1e-5))
                        @show n, b, abs.(((f  ∘ invf).(testvec) .- testvec)./testvec)
                    end
                    @test (f  ∘ invf).(testvec) ≈ testvec rtol=1e-5
                end
            end
        end

    end

    @testset "Large argument" begin

        # Some large values to test our large argument asymptotics
        large_values = [1.0e8+0im, 1.0e8+0.5im, 1.0e8-0.5im]

        @testset "H ≈ H_asymptotic" begin
            @test large_arg_hankelh1.(0, large_values) ≈ hankelh1.(0, large_values) rtol=1e-4
            @test large_arg_hankelh1.(1, large_values) ≈ hankelh1.(1, large_values) rtol=1e-4
            @test large_arg_hankelh1.(2, large_values) ≈ hankelh1.(2, large_values) rtol=1e-4
            @test large_arg_hankelh1.(3, large_values) ≈ hankelh1.(3, large_values) rtol=1e-4
            @test large_arg_hankelh1.(4, large_values) ≈ hankelh1.(4, large_values) rtol=1e-4
            @test large_arg_hankelh1.(8, large_values) ≈ hankelh1.(8, large_values) rtol=1e-4
        end

        # Test whether f f⁻¹ = identity for different Hankel order (n) and solution branches (b)
        # We apply f⁻¹ first because action of f elliminates effect of branch
        @testset "f f⁻¹(z) = z" begin
            for n in [0,1,2,4,8,11]
                for b in [0,1,2,3]
                    f   = z->   large_arg_hankelh1(n, z)
                    invf = z->inv_large_arg_hankelh1(n, z, b)
                    # f is large argument asymptotic, so f(large_values) should be in domain of invf
                    testvec = f.(large_values)
                    diffvec = abs.((f  ∘ invf).(testvec) .- testvec)
                    if !all(isapprox.((f  ∘ invf).(testvec), testvec; rtol=1e-5))
                        @show n, b, testvec, invf.(testvec), (f  ∘ invf).(testvec)
                    end
                    @test (f  ∘ invf).(testvec) ≈ testvec rtol=1e-5
                end
            end
        end

    end

end

@testset "Inverse Normalised Hankel function" begin

    @testset "Sorted Vec optimisation" begin

        hns = 1.0:-0.01:0.01
        for z₀ in [1.0, 2.0, 4.0]
            for ν in [0, 1, 4]
                # Test that the optimised and unoptimise dversion give the same result
                @test invhankelh1n.(ν, z₀, hns) ≈ invhankelh1n_sortedvec(ν, z₀, hns)
                # Test that optimised version is faster
                b1 = @benchmark invhankelh1n.($ν, $z₀, $hns) samples=10 evals=1
                b2 = @benchmark invhankelh1n_sortedvec($ν, $z₀, $hns) samples=10 evals=1
                @test minimum(b2.times) <= minimum(b1.times)
            end
        end
    end

    @testset "Householder corrector orders" begin

        # Collect the total number of predictor and corrector steps so that we can assess performance
        total_predictor_steps = 0
        total_corrector_steps1 = 0
        total_corrector_steps2 = 0
        total_corrector_steps3 = 0

        for ν in [0, 1, 4]
            for z₀ in [1.0, 2.0, 4.0]
                for ξ_target in [0.9, 0.34, 0.01]
                    iters1 = []; iters2 = []; iters3 =[];
                    ans1 = invhankelh1n_adaptive_solve(ν, z₀, ξ_target; householder_order=1, iter_vec=iters1)
                    ans2 = invhankelh1n_adaptive_solve(ν, z₀, ξ_target; householder_order=2, iter_vec=iters2)
                    ans3 = invhankelh1n_adaptive_solve(ν, z₀, ξ_target; householder_order=3, iter_vec=iters3)

                    # All methods should give similar answers
                    @test all(ans1 .≈ ans2)
                    @test all(ans2 .≈ ans3)
                    @test all(ans3 .≈ ans1)

                    # Higher order methods should require no more steps
                    @test all(iters2 .<= iters1)
                    @test all(iters3 .<= iters2)

                    # All methods should use the same predictor
                    @test length(iters1) == length(iters2) && length(iters2) == length(iters3)

                    total_predictor_steps += length(iters1)
                    total_corrector_steps1 += sum(iters1)
                    total_corrector_steps2 += sum(iters2)
                    total_corrector_steps3 += sum(iters3)

                    # Unimplemented/invalid householder corrector order
                    @test_throws Exception invhankelh1n_adaptive_solve(ν, z₀, ξ_target; householder_order=0)
                    @test_throws Exception invhankelh1n_adaptive_solve(ν, z₀, ξ_target; householder_order=4)
                end
            end
        end
        # Ratios of corrector to predictor steps for each method order
        cp_ratio1 = total_corrector_steps1/total_predictor_steps
        cp_ratio2 = total_corrector_steps2/total_predictor_steps
        cp_ratio3 = total_corrector_steps3/total_predictor_steps

        # Higher order methods should require less steps
        @test cp_ratio3 < cp_ratio2 < cp_ratio1

        # Empirically derived  differences in steps, these may change if the examples change
        @test cp_ratio2 <= cp_ratio1*0.7
        @test cp_ratio3 <= cp_ratio2*0.85
    end

    @testset "Too many steps kill switch" begin
        @test typeof(invhankelh1n_adaptive_solve(0, 1.0, 0; silent_failure=true)) <: Tuple
        @test_throws Exception invhankelh1n_adaptive_solve(0, 1.0, 0; silent_failure=false)
    end

end

@testset "Derivatives of Hankel functions   " begin

    for ν in [0, 1, 4, 10]
        for z in [1.0+0.0im, -1.0+0.1im, 0.0+2.0im, 0.0-2.0im, 100.0+5.0im]

            findiff(f, δ=1e-8) = (f(ν, z+δ) .- f(ν, z))./δ

            # Cache values
            h = hankelh1(ν, z)
            hm1 = hankelh1(ν-1, z)

             hp1               = diffhankelh1( ν, z, h, hm1)
            (hp2, hpp2)        = diff2hankelh1(ν, z, h, hm1)
            (hp3, hpp3, hppp3) = diff3hankelh1(ν, z, h, hm1)

            # Test they are consistent
            @test hp1 ≈ hp2
            @test hp2 ≈ hp3
            @test hp3 ≈ hp1
            @test hpp2 ≈ hpp3

            # Test optional arguments don't matter
            @test diffhankelh1(ν, z, h, hm1) == diffhankelh1(ν, z)
            @test diff2hankelh1(ν, z, h, hm1) == diff2hankelh1(ν, z)
            @test diff3hankelh1(ν, z, h, hm1) == diff3hankelh1(ν, z)

            # Test the derivatives are similar to naive finite difference
            @test diffhankelh1(ν,z)[1]  ≈ findiff(hankelh1)         rtol=1e-6
            @test diff2hankelh1(ν,z)[2] ≈ findiff(diffhankelh1)[1]  rtol=1e-6
            @test diff3hankelh1(ν,z)[3] ≈ findiff(diff2hankelh1)[2] rtol=1e-6
        end
    end

end
