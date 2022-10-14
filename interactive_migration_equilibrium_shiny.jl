
using Pkg
Pkg.add(["Plots", "Mux", "Interact"])

using Plots
using Mux
using Interact

#WT into GD equilbrium heatmap
mp = @manipulate throttle= 0.1 for  c_input = Widgets.slider((0:0.01:1.0), value = 0.01, label = "Conversion Efficiency"),
                                    F_input = Widgets.slider((0:0.01:1.0), value = 0.0, label = "Inbreeding")





    range1 = -5:0.1:2
    range2 = 0.001:0.01:0.99
    upper = fill(1.0im, length(range1), length(range2))
    lower = fill(1.0im, length(range1), length(range2))

    for pm_value in 1:length(range1)
            for s_value in 1:length(range2)
                    for h_value in 1
                            for c_value in c_input
                                for F_value in F_input

                                    pos1 = pm_value
                                    pos2 = s_value

                                    pm = 10^range1[pm_value]
                                    s = range2[s_value]
                                    h = h_value
                                    c = c_value
                                    F = F_value


                                    upper[pos1, pos2] = (2*c*F*h*pm*s+sqrt(Complex((-2*c*F*h*pm*s-3*c*F*h*s+2*c*F*pm*s+4*c*F*s+c*(-F)+2*c*h*pm*s+3*c*h*s-2*c*pm*s-4*c*s+c+2*F*h*pm*s+3*F*h*s-F*pm*s-2*F*s-2*h*pm*s-3*h*s+s)^2-4*(c*F*h*s-2*c*F*s+c*F-c*h*s+2*c*s-c-F*h*s+F*s+h*s+pm)*(2*c*F*h*pm*s+2*c*F*h*s-2*c*F*pm*s-2*c*F*s-2*c*h*pm*s-2*c*h*s+2*c*pm*s+2*c*s-2*F*h*pm*s-2*F*h*s+F*pm*s+F*s+2*h*pm*s+2*h*s-pm*s-s)))+3*c*F*h*s-2*c*F*pm*s-4*c*F*s+c*F-2*c*h*pm*s-3*c*h*s+2*c*pm*s+4*c*s-c-2*F*h*pm*s-3*F*h*s+F*pm*s+2*F*s+2*h*pm*s+3*h*s-s)/(2*(F-1)*(pm+1)*s*(2*c*h-2*c-2*h+1))
                                    lower[pos1, pos2] = (2*c*F*h*pm*s-sqrt(Complex((-2*c*F*h*pm*s-3*c*F*h*s+2*c*F*pm*s+4*c*F*s+c*(-F)+2*c*h*pm*s+3*c*h*s-2*c*pm*s-4*c*s+c+2*F*h*pm*s+3*F*h*s-F*pm*s-2*F*s-2*h*pm*s-3*h*s+s)^2-4*(c*F*h*s-2*c*F*s+c*F-c*h*s+2*c*s-c-F*h*s+F*s+h*s+pm)*(2*c*F*h*pm*s+2*c*F*h*s-2*c*F*pm*s-2*c*F*s-2*c*h*pm*s-2*c*h*s+2*c*pm*s+2*c*s-2*F*h*pm*s-2*F*h*s+F*pm*s+F*s+2*h*pm*s+2*h*s-pm*s-s)))+3*c*F*h*s-2*c*F*pm*s-4*c*F*s+c*F-2*c*h*pm*s-3*c*h*s+2*c*pm*s+4*c*s-c-2*F*h*pm*s-3*F*h*s+F*pm*s+2*F*s+2*h*pm*s+3*h*s-s)/(2*(F-1)*(pm+1)*s*(2*c*h-2*c-2*h+1))
                                end

                            end
                    end
            end
    end

    #upper is +complex, which ends up being less?

    realupper = real.(upper)
    realupper[imag.(upper) .!= 0] .= NaN

    reallower = real.(lower)
    reallower[imag.(lower) .!= 0] .= NaN


    realupper[realupper .< 0] .= 0
    realupper[realupper .> 1] .= NaN
    reallower[reallower .< 0] .= 0
    reallower[reallower .> 1] .= NaN

    #finding points in both equations that give desired outcome


    h1 = Plots.heatmap(range2, range1, realupper, xlab = "s", ylab = "m 10^", title = "Stable")

    h2 = Plots.heatmap(range2, range1, reallower, xlab = "s", ylab = "m 10^", title = "Unstable")

    Plots.plot(h1,h2, plot_title = "WT into GD population")

                                    
end

@layout! mp vbox(hbox(:c_input, :F_input), hbox(observe(mp)))


ui = dom"div"(mp);
WebIO.webio_serve(page("/", req -> ui), 8002);


