using Plots
using InverseHankelFunction
using LaTeXStrings

include("complex_plot_utils.jl")

function plot_hankelh1n_real_arg(ν, z₀; save=false, build_dir=".")

    zs = vcat(10.0 .^(-2:0.05:0.05), 1.0:0.1:30.0)
    hbars = hankelh1n(ν, z₀, hbars)

    plot(zs, hcat(real(hbars),imag(hbars)), label=["Real" "Imaginary"], ylab="\$\\mathrm{invhankelh1n}($ν, $z₀, \\hat{h})\$", xlab=L"\hat{h}")
    hline!([1.0], label="", linestyle=:dash)
    annotate!(1, (real(z₀) + imag(z₀))/2, text("\$\\mathrm{hankelh1n}($ν, $z₀, $z₀) = 1.0\$", :bottom))

    if save
        savefig(joinpath(build_dir,"hankelh1n_nu_$(ν)_z_0_$(z₀)_real_arg.svg"))
    end

    plot!()
end

function plot_hankelh1n_complex_arg(ν, z₀; lim=5, n_samples=200, save=false, build_dir=".")
    plot_angle_heatmap(z->hankelh1n(ν,z₀,z); lim=lim, n_samples=n_samples)

    hbars = vcat(10.0 .^(-5:0.1:0.05), 1.0:0.1:30.0)
    zs = invhankelh1n(ν, z₀, hbars)
    plot!(real.(zs), imag.(zs), label="", color=:white)
    plot!(colorbar_title="\$\\mathrm{arg}(\\mathrm{H}_{$ν}^{(1)}(z)/\\mathrm{H}_{$ν}^{(1)}($z₀))\$")
    scatter!([real.(z₀)], [imag.(z₀)], label="")
    if save
        plot!(size=(800,600))
        savefig(joinpath(build_dir,"hankelh1n_nu_$(ν)_z_0_$(z₀)_complex_arg.png"))
    end
    plot!()
end


function plot_invhankelh1n_real_arg(ν, z₀; save=false, build_dir=".")

    hbars = vcat(10.0 .^(-2:0.05:0.05), 1.0:0.1:30.0)
    zs = invhankelh1n(ν, z₀, hbars)

    plot(hbars, hcat(real(zs),imag(zs)), label=["Real" "Imaginary"], ylab="\$\\mathrm{invhankelh1n}($ν, $z₀, \\hat{h})\$", xlab=L"\hat{h}")
    vline!([1.0], label="", linestyle=:dash)
    annotate!(1, (real(z₀) + imag(z₀))/2, text("\$\\mathrm{invhankelh1n}($ν, $z₀, 1) = $z₀\$", :left))

    if save
        savefig(joinpath(build_dir,"invhankelh1n_nu_$(ν)_z_0_$(z₀)_real_arg.svg"))
    end

    plot!()
end

function make_hankelh1n_plots()
    plots_dir = (pwd()[end-3:end] == "docs") ? "build/plots" : "docs/build/plots"
    mkpath(plots_dir)

    plot_hankelh1n_complex_arg(0, 1; n_samples=1000, save=true, build_dir=plots_dir)
    plot_invhankelh1n_real_arg(0, 1; save=true, build_dir=plots_dir)

    plot_hankelh1n_complex_arg(3, 2; n_samples=1000, save=true, build_dir=plots_dir)
    plot_invhankelh1n_real_arg(3, 2; save=true, build_dir=plots_dir)
end
