module FourBars

# export compute_generic_solutions, save_generic_solutions, load_generic_solutions,
#     solve_instance, four_bars, FourBar, animate, configurations
export start_demo


import HomotopyContinuation
const HC = HomotopyContinuation
using LinearAlgebra
using Parameters: @unpack
using AbstractPlotting, GLMakie
using AbstractPlotting.MakieLayout
using ColorSchemes, Printf
using JLD2

# This follows the description in
#  Wampler, C. W., Morgan, A. P., and Sommese, A. J. (March 1, 1992).
#  "Complete Solution of the Nine-Point Path Synthesis Problem for Four-Bar Linkages."
#  ASME. J. Mech. Des. March 1992; 114(1): 153–159. https://doi.org/10.1115/1.2916909
struct FourBar
    x::ComplexF64
    a::ComplexF64
    y::ComplexF64
    b::ComplexF64
    P₀::ComplexF64 # coupler point
end
function Base.getproperty(F::FourBar, field::Symbol)
    if field === :A
        getfield(F, :P₀) + getfield(F, :a)
    elseif field === :B
        getfield(F, :P₀) + getfield(F, :b)
    elseif field === :C
        getfield(F, :P₀) + getfield(F, :y)
    elseif field === :D
        getfield(F, :P₀) + getfield(F, :x)
    elseif field === :u
        getfield(F, :x) - getfield(F, :a)
    elseif field === :v
        getfield(F, :y) - getfield(F, :b)
    else
        getfield(F, field)
    end
end

function FourBar(;
    A::ComplexF64,
    B::ComplexF64,
    C::ComplexF64,
    D::ComplexF64,
    P₀ = (D + C) / 2,
)
    x = D - P₀
    a = A - P₀
    y = C - P₀
    b = B - P₀
    FourBar(x, a, y, b, P₀)
end

function FourBar(F::FourBar, angles)
    FourBar(
        A = F.A,
        B = F.B,
        C = (F.B + F.v * angles.μ),
        D = (F.A + F.u * angles.λ),
        P₀ = (F.A + F.u * angles.λ - F.x * angles.θ),
    )
end

function four_bar_from_lengths(; A::ComplexF64, B::ComplexF64, BC, CD, AD)
    HC.@var a[1:2] b[1:2] c[1:2] d[1:2] l[1:5]
    F = HC.System(
        [
            sum((b - c) .^ 2) - BC^2,
            sum((c - d) .^ 2) - CD^2,
            sum((a - d) .^ 2) - AD^2,
            sum(l .* [c; d; 1]),
        ],
        [c; d],
        [l; a; b],
    )

    for i = 1:5
        res =
            HC.solve(F, target_parameters = [randn(5); real(A); imag(A); real(B); imag(B)])
        S = HC.real_solutions(res)
        if !isempty(S)
            s = S[1]
            C = complex(s[1], s[2])
            D = complex(s[3], s[4])
            return FourBar(A = A, B = B, C = C, D = D)
        end
    end
    error("Could not find a valid configuartion.")
end

Base.show(io::IO, F::FourBar) = print_fieldnames(io, F)

"""
     print_fieldnames(io::IO, obj)

 A better default printing for structs.
 """
function print_fieldnames(io::IO, obj)
    println(io, typeof(obj), ":")
    for name in fieldnames(typeof(obj))
        if getfield(obj, name) !== nothing
            println(io, " • ", name, " → ", getfield(obj, name))
        end
    end
end


is_conjugated_pair(u, v, tol) = abs(u - conj(v)) < tol
is_valid_loop_solution(μ, θ, τ, τ̂; tol::Float64 = 1e-12) =
    is_conjugated_pair(τ, τ̂, tol) && abs(abs(μ) - 1) < tol && abs(abs(θ) - 1) < tol

function loop_homotopy()
    params = HC.@var x a y b x̂ â ŷ b̂
    HC.@var λ μ θ τ τ̂
    HC.Homotopy(
        [
            (x - a) * λ - x * θ - τ + a,
            (x̂ - â) * θ - x̂ * λ - (τ̂ - â) * λ * θ,
            (y - b) * μ - y * θ - τ + b,
            (ŷ - b̂) * θ - ŷ * μ - (τ̂ - b̂) * μ * θ,
        ],
        [μ, θ, τ, τ̂],
        λ,
        collect(params),
    )
end

function trace_curve(F::FourBar; Δt = 1e-2, max_steps = 1_000, accuracy = 1e-8)
    H = loop_homotopy()
    params = let
        @unpack x, a, y, b = F
        x̂, â, ŷ, b̂ = conj.((x, a, y, b))
        [x, a, y, b, x̂, â, ŷ, b̂]
    end
    # compute once solutions for our loop homotopy at λ₀
    λ_gen = randn(ComplexF64)
    F = HC.System(H.expressions, H.variables, [H.t; H.parameters])
    generic_solutions = HC.solutions(HC.solve(F, target_parameters = [λ_gen; params]))


    # a valid solution is
    φ = φ₀ = 0.0
    μ, θ, τ, τ̂ = [1.0, 1, 0, 0im]

    angles = [(λ = cis(φ), μ = μ, θ = θ)]
    # loop, (λ, _) = loop_equations(F)
    # φ = φ₀ = angle(λ₀)
    # angles = [(cis(φ₀), x₀[1], x₀[2])]
    # μ₀, θ₀ = angle(x₀[1]), angle(x₀[2])
    tracker = HC.Tracker(HC.fix_parameters(H, params); options = (max_steps = 200,))

    # tracker = coretracker(SPHomotopy(loop, λ), [randn(ComplexF64, 4)]; accuracy = accuracy)
    HC.init!(tracker, [1.0, 1, 0, 0im], cis(φ), cis(φ + Δt))
    x = tracker.state.x
    y = copy(x)
    Δ = Δt
    for i = 2:max_steps
        if φ * (φ + Δ) < 0 # jump over zero
            φ′ = 0.0
        else
            φ′ = φ + Δ
        end
        retcode = HC.track!(tracker, y, cis(φ), cis(φ′))
        if retcode != HC.TrackerCode.success
            @warn "PathTracker failed with $retcode"
            break
        end
        μ, θ, τ, τ̂ = x
        if is_valid_loop_solution(μ, θ, τ, τ̂)
            if abs(angle(angles[end].μ) - angle(μ)) < 2Δt
                φ = φ′
                Δ = clamp(2Δ, -Δt, Δt)
                push!(angles, (λ = cis(φ), μ = μ, θ = θ))
                y .= x
                # reduce step size since μ jump is too large
            else
                Δ /= 2
                continue
            end
        else
            # jump to different branch
            branch_solutions = map(generic_solutions) do s
                HC.solution(HC.track(tracker, s, λ_gen, cis(φ)))
            end
            y .= branch_solutions[last(findmax(norm.(branch_solutions .- Ref(y))))]
            Δ = -Δ
            μ, θ = y[1], y[2]
            push!(angles, (λ = cis(φ), μ = μ, θ = θ))
        end
        if φ > π
            φ -= 2π
        elseif φ < -π
            φ += 2π
        end
        if abs(φ - φ₀) < 0.5Δt && abs(1 - y[1]) < 1e-3 && abs(1 - y[2]) < 1e-3
            pop!(angles)
            break
        end
    end
    angles
end



fourbar_to_points(F) = Point.(reim.((F.A, F.B, F.C, F.D)))


function energy(
    fourbar;
    resting_lengths::Vector{Float64},
    elasticities::Vector{Float64},
    E::ComplexF64,
    F::ComplexF64,
)
    c1, c2 = elasticities
    r1, r2 = resting_lengths

    c1 * (max(0.0, sqrt(abs2(fourbar.C - E)) - r1))^2 +
    c2 * (max(0.0, sqrt(abs2(fourbar.D - F)) - r2))^2
end

function energy_system(;
    bar_lengths::Vector,
    resting_lengths,
    elasticities,
    A::ComplexF64,
    B::ComplexF64,
)
    HC.@var p[1:2, 1:6]

    bars = [(2, 3), (3, 4), (4, 1)]
    cables = [(3, 5), (4, 6)]
    HC.@var δ[1:2] λ[1:5]

    G = [
        [
            sum((p[:, i] - p[:, j]) .^ 2) .- bar_lengths[k]^2
            for (k, (i, j)) in enumerate(bars)
        ]
        [sum((p[:, i] - p[:, j]) .^ 2) .- δ[k]^2 for (k, (i, j)) in enumerate(cables)]
    ]
    Q = sum(elasticities .* (δ .- resting_lengths) .^ 2)
    L = Q + λ'G

    X = [p[:, 3]; p[:, 4]]
    # fix constants
    L′ = HC.subs(L, p[:, 1] => [real(A), imag(A)], p[:, 2] => [real(B), imag(B)])
    Y = [p[:, 5]; p[:, 6]]
    ∇L = HC.differentiate(L′, [X; δ; λ])
    F = HC.System(∇L, variables = [X; δ; λ], parameters = Y)
end

function energy_discriminant_system(; kwargs...)
    F = energy_system(; kwargs...)
    exprs = HC.expressions(F)
    Y = F.parameters
    HC.@var v[1:HC.nvariables(F)] a[1:2] b
    normalizer = [
        0.7459631719784344,
        -0.09612081503288031,
        0.7270390746748708,
        0.18865990309925695,
        -0.07948057722799588,
        -0.4765286254244125,
        1.468490025187829,
        -1.225830034198433,
        -1.0520616528796696,
        -0.28900616789894923,
        0.2838676399853993,
    ]
    H = HC.System(
        [
            exprs
            HC.differentiate(exprs, F.variables) * v
            normalizer'v - 1
            a' * Y[3:4] + b
        ],
        variables = [Y[3:4]; F.variables; v],
        parameters = [Y[1:2]; a; b],
    )
end

function local_minimum_checker(;
    bar_lengths::Vector,
    resting_lengths,
    elasticities,
    A::ComplexF64,
    B::ComplexF64,
)
    HC.@var p[1:2, 1:6]

    bars = [(2, 3), (3, 4), (4, 1)]
    cables = [(3, 5), (4, 6)]

    HC.@var l²[1:3] r[1:2] δ[1:2] λ[1:5]

    G = HC.subs(
        [
            [
                sum((p[:, i] - p[:, j]) .^ 2) .- bar_lengths[k]^2
                for (k, (i, j)) in enumerate(bars)
            ]
            [sum((p[:, i] - p[:, j]) .^ 2) .- δ[k]^2 for (k, (i, j)) in enumerate(cables)]
        ],
        p[:, 1] => [real(A), imag(A)],
        p[:, 2] => [real(B), imag(B)],
    )

    Q = sum(elasticities .* (δ .- resting_lengths) .^ 2)
    L = Q + λ'G

    X = [p[:, 3]; p[:, 4]]
    # fix constants
    L′ = HC.subs(L, p[:, 1] => [real(A), imag(A)], p[:, 2] => [real(B), imag(B)])
    Y = [p[:, 5]; p[:, 6]]
    ∇L = HC.differentiate(L′, [X; δ])
    HL = HC.CompiledSystem(HC.System(∇L, [X; δ], [λ; Y]))
    dG = HC.CompiledSystem(HC.System(G, [X; δ], [λ; Y]))

    (s, params) -> begin
        v = s[1:6]
        q = [s[7:end]; params]
        W = HC.jacobian!(zeros(6, 6), HL, v, q)
        V = nullspace(HC.jacobian!(zeros(5, 6), dG, v, q))
        E = eigvals(V' * W * V)
        all(e -> e > 1e-8, E)
    end
end




function comp_min_energy_positions(;
    A::ComplexF64,
    B::ComplexF64,
    E::ComplexF64,
    F::ComplexF64,
    bar_lengths,
    elasticities,
    resting_lengths,
)

    energy_sys = energy_system(
        A = A,
        B = B,
        bar_lengths = bar_lengths,
        elasticities = elasticities,
        resting_lengths = resting_lengths,
    )
    gen_params = randn(ComplexF64, 4)
    gen_res = HC.solve(energy_sys; target_parameters = gen_params)


    min_checker = local_minimum_checker(
        A = A,
        B = B,
        bar_lengths = bar_lengths,
        elasticities = elasticities,
        resting_lengths = resting_lengths,
    )
    res = HC.solve(
        energy_sys,
        gen_res;
        start_parameters = gen_params,
        target_parameters = [reim(E)..., reim(F)...],
    )

    min_energy_sols = filter(
        s -> s[5] > -1e-12 && s[6] > -1e-12 && min_checker(s, [reim(E)..., reim(F)...]),
        HC.real_solutions(res),
    )
    min_energy_positions = map(min_energy_sols) do s
        FourBar(; A = A, B = B, C = s[1] + im * s[2], D = s[3] + im * s[4])
    end
end

# coordinate transformation between two rects
function transferrects(pos, rectfrom, rectto)
    fracpos = (pos .- rectfrom.origin) ./ rectfrom.widths
    fracpos .* rectto.widths .+ rectto.origin
end

function add_control_node!(ax, n; kwargs...)
    selected = Ref(false)
    plt = scatter!(ax, n; kwargs...)
    lift(events(ax.scene).mouseposition) do pos
        x, y = transferrects(pos, ax.scene.px_area[], ax.limits[])
        if AbstractPlotting.is_mouseinside(ax.scene) && selected[]
            p = Point2f0(x, y)
            n[] = p
        end
        nothing
    end
    on(events(ax.scene).mousedrag) do drag
        if ispressed(ax.scene, Mouse.left)
            if drag == Mouse.down
                plot, _idx = mouse_selection(ax.scene)
                if plot == plt
                    selected[] = true
                end
            end
        elseif drag == Mouse.up || !AbstractPlotting.is_mouseinside(ax.scene)
            selected[] = false
        end
    end
    n
end

function compute_catastrophe_set()
    H = energy_discriminant_system(
        A = A,
        B = B,
        bar_lengths = bar_lengths,
        elasticities = elasticities,
        resting_lengths = resting_lengths,
    )
    mres = HC.monodromy_solve(H, target_solutions_count = 288, max_loops_no_progress = 10)
    N = 600
    xpoints = HC.solve(
        H,
        HC.solutions(mres),
        start_parameters = HC.parameters(mres),
        target_parameters = vcat.(reim(E)..., 1, 0, range(-6.5; stop = 6.5, length = N)),
        transform_result = (r, p) -> HC.real_solutions(r),
        flatten = true,
    )
    ypoints = HC.solve(
        H,
        HC.solutions(mres),
        start_parameters = HC.parameters(mres),
        target_parameters = vcat.(reim(E)..., 0, 1, range(-6.5; stop = 6.5, length = N)),
        transform_result = (r, p) -> HC.real_solutions(r),
        flatten = true,
    )

    points = [xpoints; ypoints]
    rpoints = filter(p -> p[7] ≥ -0 && p[8] ≥ 0, points)

    (discriminant = map(p -> p[1:2], points), catastrophe_set =  map(p -> p[1:2], rpoints))
end



function start_demo(;
    A = -1.0 + 0im,
    B = 1.0 + 0im,
    E = 4.0 + 3im,
    F = -6.0 + 2im,
    bar_lengths = [3.0, 2.0, 1.5],
    elasticities = [1.0, 2.0],
    resting_lengths = [0.1, 0.1],
    # discriminant = nothing,
    catastrophe_set = nothing,
)

    fourbar_data = jldopen(joinpath(@__DIR__, "../data/fourbar_catastrophe.jld2"))
    # if discriminant === nothing
    #     discriminant = fourbar_data["discriminant"]
    # end
    if catastrophe_set === nothing
        catastrophe_set = fourbar_data["catastrophe_set"]
    end

    scene, layout = layoutscene(30, resolution = ((1600, 2000)), font = "KaTeX_Main")
    scene.fontsize = 40f0
    ax = (layout[1:4, 1] = LAxis(scene))
    ax.ylabel = "y"
    ax.xlabel = "x"

    energy_ax = (layout[5, 1] = LAxis(
        scene,
        title = "Potential energy Q along coupler curve",
        # titlefont = "SanFrancisco",
    ))
    ax.xgridvisible = ax.ygridvisible = false
    energy_ax.xgridvisible = energy_ax.ygridvisible = false
    energy_ax.ylabel = "Q"
    energy_ax.ytickformat = xs -> [@sprintf("%03.0f", x) for x in xs]
    # energy_ax.ytickformat = xs -> [@sprintf("%.0f", x) for x in xs]
    hidexdecorations!(energy_ax)


    min_energy_positions = comp_min_energy_positions(;
        A = A,
        B = B,
        E = E,
        F = F,
        bar_lengths = bar_lengths,
        elasticities = elasticities,
        resting_lengths = resting_lengths,
    )

    fourbar = min_energy_positions[1]
    angles = trace_curve(fourbar; Δt = 1e-2, max_steps = 10_000)
    loop = FourBar.(Ref(fourbar), angles)


    loop_index = Node(1)
    NABCD = lift(k -> fourbar_to_points(loop[k]), loop_index)
    NF = Node(Point(reim(F)))


    N_energy = lift(NF) do f
        energy.(
            loop;
            E = E,
            F = complex(f[1], f[2]),
            elasticities = elasticities,
            resting_lengths = resting_lengths,
        )
    end

    on(N_energy) do energy
        n = length(energy)
        k = loop_index[]
        # let's go to the right
        # println(mod1(k - 1, n), " ", k, " ", mod1(k + 1, n))
        # println(energy[mod1(k - 1, n)], " ", energy[k], " ", energy[mod1(k + 1, n)])
        while energy[mod1(k + 1, n)] < energy[k] #* (1 + 1e-14)
            k = mod1(k + 1, n)
        end
        # @show energy[mod1(k + 1, n)] < energy[k]
        # if we didn't do a step, go to the left
        if k == loop_index[]
            while energy[mod1(k - 1, n)] < energy[k]# * (1 + 1e-14)
                k = mod1(k - 1, n)
            end
        end
        loop_index[] = k
    end

    lines!(
        ax,
        map(fb -> Point(reim(fb.P₀)), loop),
        linestyle = :dash,
        linewidth = 4,
        color = :gray,
    )


    scatter!(
        ax,
        Point{2,Float64}.(catastrophe_set),
        markersize = 4px,
        color = :LIGHTCORAL,
        strokewidth = 0,
    )
    # scene = Scene(limits = limits, resolution = (1200, 1200), scale_plot = false)
    linesegments!(ax, lift(NABCD) do (A, B, C, D)
        [B => C, C => D, D => A]
    end, linewidth = 4, color = :black)
    linesegments!(
        ax,
        lift(NABCD, NF) do (A, B, C, D), F
            [D => F]
        end,
        linewidth = 4,
        color = :MEDIUMAQUAMARINE,
    )
    linesegments!(ax, lift(NABCD) do (A, B, C, D)
        [C => Point(reim(E))]
    end, linewidth = 4, color = :MEDIUMAQUAMARINE)

    add_control_node!(
        ax,
        NF;
        markersize = 20px,
        marker = :cross,
        color = :black,
    )
    scatter!(ax, lift(x -> collect(x)[3:4], NABCD), color = :black)
    scatter!(ax, lift(x -> collect(x)[1:2], NABCD), color = :black, marker = :rect)
    scatter!(ax, Point(reim(E)), marker = :rect, color = :black)


    NABCD[] = fourbar_to_points(first(loop))

    xlims!(ax, [-6.5, 6.5]) # as vector
    ylims!(ax, [-6.5, 6.5]) # as vector
    ax.aspect = AxisAspect(1)

    loop_cum_length = cumsum(map(enumerate(loop)) do (k, fb)
        abs(fb.P₀ - loop[max(k - 1, 1)].P₀)
    end)

    lines!(energy_ax, loop_cum_length, N_energy, linewidth = 4)
    scatter!(
        energy_ax,
        lift((i, E) -> Point(loop_cum_length[i], E[i]), loop_index, N_energy),
        color = :MEDIUMAQUAMARINE,
        markersize = 24px,
        marker = :diamond,
        strokewidth = 0,
    )
    autolimits!(energy_ax)
    xlims!(energy_ax, [loop_cum_length[1], loop_cum_length[end]])
    on(N_energy) do _
        autolimits!(energy_ax)
        xlims!(energy_ax, [loop_cum_length[1], loop_cum_length[end]])
    end
    display(scene)
    scene
end
end #module
