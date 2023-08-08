module App
using GenieFramework, StippleLatex, DifferentialEquations, ModelingToolkit
@genietools

function define_ODE()
    @parameters t σ ρ β
    @variables x(t) y(t) z(t)
    D = Differential(t)

    eqs = [D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z]
    @named sys = ODESystem(eqs)
    sys = structural_simplify(sys)
    u0 = [x => 1.0,
        y => 0.0,
        z => 0.0]
    p = [σ => 10.0,
        ρ => 28.0,
        β => 8 / 3]

    tspan = (0.0, 100.0)
    ODEProblem(sys, u0, tspan, p, jac=true)
end

prob = define_ODE()


@handlers begin
    @private integrator = DifferentialEquations.init(prob, Tsit5())
    @in σ = 10
    @in ρ = 28.0
    @in β = 8 / 3
    @out t::Float32 = 0.0
    @in t_step = 0.03
    @in t_end = 10
    @in start = false
    @out solplot = PlotData()
    @out layout = PlotLayout(
        xaxis=[PlotLayoutAxis(xy="x", title="x")],
        yaxis=[PlotLayoutAxis(xy="y", title="y")])
    @private u_x = []
    @private u_y = []
    @private running = false
    @onchange start begin
        DifferentialEquations.reinit!(integrator)
        u_x = []
        u_y = []
        t = 0.0
        if running == false
            running = true
            @async begin
                while t <= t_end
                    sleep(t_step)
                    solplot = PlotData(x=u_x, y=u_y, z=u_x, plot=StipplePlotly.Charts.PLOT_TYPE_LINE)
                    integrator.p[1] = σ
                    integrator.p[2] = ρ
                    integrator.p[3] = β
                    step!(integrator, t_step, true)
                    l = length(integrator.sol.u)
                    append!(u_x, [u[1] for u in integrator.sol.u[l:end]][:])
                    append!(u_y, [u[2] for u in integrator.sol.u[l:end]][:])
                    t = integrator.sol.t[end]
                end
                running = false
            end
        end
    end
end


meta = Dict("og:title" => "Lorenz Chaotic Attractor", "og:description" => "Real-time simulation of a dynamic system with constant UI refresh.", "og:image" => "/preview.jpg")
layout = DEFAULT_LAYOUT(meta=meta)
@page("/", "app.jl.html", layout)
end
