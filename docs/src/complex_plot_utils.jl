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