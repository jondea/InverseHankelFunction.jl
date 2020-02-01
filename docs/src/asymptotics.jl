using Plots
splitreim(v) = [real(v) imag(v)]

function plot_hankel_with_approx(ν, dir)
    x = 0.01:0.01:10.0

    plot(grid=false, xlab="x", ylab="f", ylims=(-3,1), legend=:bottomright, title="hankelh1 (order $ν) and its small/large argument asymptotics")

    plot!(x, splitreim(hankelh1.(ν, x)), color=[:orange :blue], linestyle=:solid, label=["Real" "Imag"])

    x_small = 0.01:0.01:2.0
    plot!(x_small, splitreim(small_arg_hankelh1.(ν, x_small)), color=[:orange :blue], linestyle=:dash, label=["Real small arg" "Imag small arg"])

    x_large = 1.0:0.01:10.0
    plot!(x_large, splitreim(large_arg_hankelh1.(ν, x_large)), color=[:orange :blue], linestyle=:dot, label=["Real large arg" "Imag large arg"])

    savefig("$dir/hankel_with_approx_nu_$ν.svg")
end

function make_asymptotic_plots()

    plots_dir = (pwd()[end-3:end] == "docs") ? "build/plots" : "docs/build/plots"
    mkpath(plots_dir)

    plot_hankel_with_approx(0, plots_dir)
    plot_hankel_with_approx(1, plots_dir)
    plot_hankel_with_approx(2, plots_dir)

end
