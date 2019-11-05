using InverseHankelFunction
import Test: @test, @testset
import SpecialFunctions: hankelh1

include("utils.jl")

@testset "Hankel function asymptotics" begin
    
    @testset "Small argument" begin

        # Some small values to test our small argument asymptotics
        small_values = 1e-4 .*[1+0im, 1im, -1+0im, -1im]

        @testset "H ≈ H_asymptotic" begin
            @test all(isapprox.(small_arg_hankelh1.(0, small_values), hankelh1.(0, small_values); rtol=0.1))
            @test all(isapprox.(small_arg_hankelh1.(1, small_values), hankelh1.(1, small_values); rtol=1e-5))
            @test all(isapprox.(small_arg_hankelh1.(2, small_values), hankelh1.(2, small_values); rtol=1e-5))
            @test all(isapprox.(small_arg_hankelh1.(4, small_values), hankelh1.(4, small_values); rtol=1e-5))
            @test all(isapprox.(small_arg_hankelh1.(8, small_values), hankelh1.(8, small_values); rtol=1e-5))
            @test all(isapprox.(small_arg_hankelh1.(11, small_values), hankelh1.(11, small_values); rtol=1e-5))
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
                    @test all(isapprox.((f  ∘ invf).(testvec), testvec; rtol=1e-5))
                end
            end
        end
        
    end

    @testset "Large argument" begin
    
        # Some large values to test our large argument asymptotics
        large_values = [1.0e8+0im, 1.0e8+0.5im, 1.0e8-0.5im]

        @testset "H ≈ H_asymptotic" begin
            @test all(isapprox.(large_arg_hankelh1.(0, large_values), hankelh1.(0, large_values); rtol=1e-4))
            @test all(isapprox.(large_arg_hankelh1.(1, large_values), hankelh1.(1, large_values); rtol=1e-4))
            @test all(isapprox.(large_arg_hankelh1.(2, large_values), hankelh1.(2, large_values); rtol=1e-4))
            @test all(isapprox.(large_arg_hankelh1.(3, large_values), hankelh1.(3, large_values); rtol=1e-4))
            @test all(isapprox.(large_arg_hankelh1.(4, large_values), hankelh1.(4, large_values); rtol=1e-4))
            @test all(isapprox.(large_arg_hankelh1.(8, large_values), hankelh1.(8, large_values); rtol=1e-4))
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
                    @test all(isapprox.((f  ∘ invf).(testvec), testvec; rtol=1e-5))
                end
            end
        end
    
    end

end

