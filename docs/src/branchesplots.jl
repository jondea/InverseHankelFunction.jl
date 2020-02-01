
using Plots
using InverseHankelFunction
using LaTeXStrings
using SpecialFunctions

function plot_branches(ν)

    # plot(legend=false, grid=false, aspectratio=1.0)
    hs = (10.0:-0.1:0.1) .+ 0.0im
    hs = ((0.01:0.01:10.0) .+ 0.0im).^-1
    for b in -ν:(ν-1)
        zs = inv_small_arg_hankelh1.(ν, hs, b)
        plot!(real.(zs), imag.(zs), linestyle=:solid, color=:cyan, label="")
    end

    for b in -ν:(ν-1)
        # initial guess chooses the branch
        z₀ = inv_small_arg_hankelh1(ν, first(hs), b)
        h₀ = first(hs)

        hbar = one(h₀)
        z = z₀

        zs = similar(hs)
        for i in eachindex(zs)
            zs[i] = invnormalisedhankel_adaptive_solve(ν, z₀, hs[i]/h₀, z, hbar; silent_failure=true)
            hbar = hs[i]/h₀
            z = zs[i]
        end
        plot!(real.(zs), imag.(zs), linestyle=:dash, color=:white, label="")
    end

    plot!()
end

function plot_q(ν, α;n_samples=200)

    function objective(z)
        q = hankelh1(ν, z)/α
        nu_zero_deriv = 1.0 - (real(q)*real(q) + imag(q)*imag(q))/real(q)
        return nu_zero_deriv
        return abs(q/(1.0-nu_zero_deriv) - 1.0)
    end

    lims = (-5,5)
    xlims = lims
    ylims = lims
    x = range(xlims[1], stop=xlims[2], length=n_samples)
    y = range(ylims[1], stop=ylims[2], length=n_samples)
    heatmap(x, y, (x,y)->log10(abs((objective(x+im*y)))))

    plot!(xlims=xlims, ylims=ylims, aspectratio=1.0, xlab=L"\mathrm{Real}(z)", ylab=L"\mathrm{Imag}(z)", colorbar_title=L"\log_{10}(q)")

end

function plot_angle_heatmap(f;lim=5, n_samples=200)
    x = range(-lim, stop=lim, length=n_samples)
    y = range(-lim, stop=lim, length=n_samples)
    heatmap(x, y, (x,y)->angle(f(x+im*y)), color=:colorwheel, xlims=(-lim,lim), ylims=(-lim,lim),
           clims=(-1π,1π), xlab=L"\mathrm{Re}(z)", ylab=L"\mathrm{Im}(z)", aspectratio=1.0)
end

function plot_angle_contour(f;lim=5, n_samples=200)
    x = range(-lim, stop=lim, length=n_samples)
    y = range(-lim, stop=lim, length=n_samples)
    contour(x, y, (x,y)->abs(angle(f(x+im*y))), color=:colorwheel, clims=(-1π,1π), xlab=L"\mathrm{Re}(z)", ylab=L"\mathrm{Im}(z)", aspectratio=1.0)
end

function plot_angle_and_abs_contours!(f;lim=5, n_samples=200)
    x = range(-lim, stop=lim, length=n_samples)
    y = range(-lim, stop=lim, length=n_samples)
    contour!(x, y, (x,y)->angle(f(x+im*y)), color=:white)
    contour!(x, y, (x,y)->log(abs(f(x+im*y))), color=:black)
end

function plot_angle_and_contours(ν;lim=5, n_samples=200, save=false, build_dir=".")
    plot_angle_heatmap(z->hankelh1(ν,z); lim=lim, n_samples=n_samples)
    plot_angle_and_abs_contours!(z->hankelh1(ν,z); lim=lim, n_samples=n_samples)
    plot!(colorbar_title="\$\\mathrm{H}_{$ν}^{(1)}(z)\$")
    if save
        plot!(size=(800,600))
        savefig(joinpath(build_dir,"hankelh1_nu_$(ν)_angle_and_contours.png"))
    end
    plot!()
end

function make_branches_plots()
    plots_dir = (pwd()[end-3:end] == "docs") ? "build/plots" : "docs/build/plots"
    mkpath(plots_dir)

    plot_angle_and_contours(0; n_samples=1000, save=true, build_dir=plots_dir)
    plot_angle_and_contours(1; n_samples=1000, save=true, build_dir=plots_dir)
    plot_angle_and_contours(2; n_samples=1000, save=true, build_dir=plots_dir)
    plot_angle_and_contours(3; n_samples=1000, save=true, build_dir=plots_dir)
    plot_angle_and_contours(8; n_samples=1000, save=true, build_dir=plots_dir)
end
