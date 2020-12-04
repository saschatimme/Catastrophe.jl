module Zeeman

export start_demo

using AbstractPlotting, GLMakie
using AbstractPlotting.MakieLayout
using HomotopyContinuation
using LinearAlgebra
using Printf

# coordinate transformation between two rects
function transferrects(pos, rectfrom, rectto)
  fracpos = (pos .- rectfrom.origin) ./ rectfrom.widths
  fracpos .* rectto.widths .+ rectto.origin
end

add_control_node!(ax, init_coords) = add_control_node!(identity, ax, init_coords)
function add_control_node!(cb, ax, n)
  selected = Ref(false)
  plt = scatter!(ax, n, markersize = 32px, marker=:cross, color=:black)
  lift(events(ax.scene).mouseposition) do pos
    x, y = transferrects(pos, ax.scene.px_area[], ax.limits[])
    if AbstractPlotting.is_mouseinside(ax.scene) && selected[]
      p = Point2f0(x, y)
      n[] = p
      cb(p)
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

function start_demo()
  F = let
    @var p[1:2, 1:2] x[1:2] c[1:2] rc r[1:2] d[1:2] λ[1:3]
    p₁, p₂ = p[:, 1], p[:, 2]

    Q̃ = sum((d - r) .^ 2) / 2
    G = [
      (sum((x - c) .^ 2) - rc^2) / 2,
      (sum((p₁ - x) .^ 2) - d[1]^2) / 2,
      (sum((p₂ - x) .^ 2) - d[2]^2) / 2,
    ]

    L = Q̃ + G' * λ
    ∇L = differentiate(L, [x; d; λ])
    F = System(∇L; variables = [x; d; λ], parameters = [p₁; p₂; r; c; rc])
  end

  params₀ = randn(ComplexF64, nparameters(F))
  res₀ = solve(F; target_parameters = params₀)
  S₀ = solutions(res₀)
  H = ParameterHomotopy(F, params₀, params₀)
  solver, _ = solver_startsolutions(H, S₀)

  is_local_minimum = let
    @var p[1:2, 1:2] x[1:2] c[1:2] rc r[1:2] d[1:2] λ[1:3]
    p₁, p₂ = p[:, 1], p[:, 2]

    Q̃ = sum((d - r) .^ 2) / 2
    G = [
      (sum((x - c) .^ 2) - rc^2) / 2,
      (sum((p₁ - x) .^ 2) - d[1]^2) / 2,
      (sum((p₂ - x) .^ 2) - d[2]^2) / 2,
    ]

    L = Q̃ + G' * λ
    ∇L = differentiate(L, [x; d])

    HL = InterpretedSystem(System(∇L, [x; d], [λ; p₁; p₂; r; c; rc]))
    dG = InterpretedSystem(System(G, [x; d], [λ; p₁; p₂; r; c; rc]))

    (s, params) -> begin
      v = s[1:4]
      q = [s[5:end]; params]
      W = zeros(4, 4)
      jacobian!(W, HL, v, q)
      J = zeros(3, 4)
      jacobian!(J, dG, v, q)
      V = nullspace(J)
      all(e -> e ≥ 1e-14, eigvals(V' * W * V))
    end
  end

  comp_next_x = let
    params₀ = randn(ComplexF64, nparameters(F))
    res₀ = solve(F; target_parameters = params₀)
    S₀ = solutions(res₀)
    H = ParameterHomotopy(F, params₀, params₀)
    solver, _ = solver_startsolutions(H, S₀)
    params = [0.0, 0.0, 2, 1 / 4, 1 / 4, 1 / 4, 1, 0, 0.5]

    (curr_x, p₁) -> begin
      best_d = Inf
      next_x = nothing
      params[1:2] .= p₁
      target_parameters!(solver, params)
      R = solve(solver, S₀, threading = false)
      # has at most two critical points
      other_critical = nothing
      for r in R
        if is_real(r)
          s = solution(r)
          if real(s[3]) ≥ 0 && real(s[4]) ≥ 0 && is_local_minimum(real.(s), params)
            next_x′ = Point2f0(real(s[1]), real(s[2]))
            d = norm(curr_x - next_x′)
            if d < best_d
              other_critical = next_x
              next_x = next_x′
              best_d = d
            else
              other_critical = next_x′
            end
          end
        end
      end
      next_x, other_critical
    end
  end

  P₁ = let
    @var p[1:2, 1:2] x[1:2] c[1:2] rc r[1:2] d[1:2] λ[1:3]
    p₁, p₂ = p[:, 1], p[:, 2]

    Q̃ = sum((d - r) .^ 2) / 2
    G = [
      (sum((x - c) .^ 2) - rc^2) / 2,
      (sum((p₁ - x) .^ 2) - d[1]^2) / 2,
      (sum((p₂ - x) .^ 2) - d[2]^2) / 2,
    ]

    L = Q̃ + G' * λ
    ∇L = differentiate(L, [x; d; λ])
    @var v[1:7]
    α = [
          0.7459631719784344,
          -0.09612081503288031,
          0.7270390746748708,
          0.18865990309925695,
          -0.07948057722799588,
          -0.4765286254244125,
          1.468490025187829,
      ]
    J_L_v = [differentiate(∇L, [x; d; λ]) * v; α'v - 1]
    @var a[1:2] b
    L₁ = a' * p₁ + b
    P₁ = System(
      [∇L; J_L_v; L₁];
      variables = [p₁; x; d; λ; v],
      parameters = [p₂; r; c; rc; a; b],
    )
  end
  monodromy_res = monodromy_solve(P₁, target_solutions_count = 72)
  rand_lin_space = let p₂ = [2, 0.25], r = [0.25, 0.25], c = [1, 0], rc = [0.5]
    fixed = [p₂; r; c; rc]
    () -> [fixed; randn(3)]
  end
  N = 400
  alg_catastrophe_points = solve(
    P₁,
    solutions(monodromy_res),
    start_parameters = parameters(monodromy_res),
    target_parameters = [rand_lin_space() for i = 1:N],
    transform_result = (r, p) -> real_solutions(r),
    flatten = true,
  )
  filter!(p -> p[5] ≥ 0 && p[6] ≥ 0, alg_catastrophe_points);
  catastrophe_points = map(p -> Point2f0(p[1], p[2]), alg_catastrophe_points)


  let
    scene, layout = layoutscene(30, resolution = (1300, 1470))
    ax = (layout[1:3, 1] = LAxis(scene, title = "Control plane"))
    ax.xzoomlock = true
    ax.yzoomlock = true
    layout[1:3, 2] = buttongrid = GridLayout(tellwidth = true)
    energy_ax =
      (layout[4, 1:2] = LAxis(scene, title = "Energy landscape", xlabel = "θ", ylabel = "Q"))

    supertitle =
      layout[0, 1] = LText(
        scene,
        "Zeeman's Catastrophe Machine with HomotopyContinuation.jl",
        textsize = 30,
        font = "Noto Sans Bold",
        color = (:black, 0.8),
      )

    p₁ = Node(Point2f0(0, 0))

    p₂ = Point2f0(2, 1 / 4)
    c = Point2f0(1, 0)

    next_x, other_critical = comp_next_x(Point2f0(0), p₁[])
    x = Node(next_x)
    if isnothing(other_critical)
      shadow_x = Node(Point2f0(0))
    else
      shadow_x = Node(other_critical)
    end

    # bars
    linesegments!(ax, @lift [c, $x]; linewidth = 6.0)
    # elastic bands
    linesegments!(
      ax,
      @lift [$p₁, $x];
      linewidth = 4,
          color = :MEDIUMAQUAMARINE,
    )
    linesegments!(ax, @lift [p₂, $x]; linewidth = 4,
        color = :MEDIUMAQUAMARINE,)

    # catastrophe curve
    catastrophe_curve = scatter!(ax, catastrophe_points, markersize = 4px, color = :indianred, strokewidth=0)

    # constants
    scatter!(ax, [p₂, c]; markersize = 32px, marker = :rect, color=:black)
    # internal
    scatter!(ax, x; color = :black, marker = :diamond, markersize = 32px)
    add_control_node!(ax, p₁)
    # other_critical
    shadow_plt = scatter!(
      ax,
      shadow_x;
      color = :lightgrey,
      transparency = true,
      alpha = 0.1,
      marker = :diamond,
      markersize = 32px,
    )

    shadow_visible = Node(!isnothing(other_critical))
    shadow_plt.visible[] = shadow_visible[]
    map!(x, p₁) do p
      next_x, other_critical = comp_next_x(x[], p)
      shadow_visible[] = !isnothing(other_critical)
      if !isnothing(other_critical)
        shadow_x[] = other_critical
      end
      next_x
    end

    xlims!(ax, [-0.75, 2.5]) # as vector
    ylims!(ax, [-1.75, 1.5]) # as vector
    ax.aspect = AxisAspect(1)

    ## energy
    θs = -π:0.001:π
    xs = 0.5 .* Point2f0.(reverse.(sincos.(θs))) .+ Ref(c)

    Q(x, p1) = 0.5 * ((norm(p1 - x) - 0.25)^2 + (norm(p₂ - x) - 0.25)^2)

    energy_landscape = lift(p1 -> map(x -> Q(x, p1), xs), p₁)
    lines!(energy_ax, θs, energy_landscape; color = :black, linewidth = 4)
    scatter!(
      energy_ax,
      lift(x, p₁) do x, p₁
        Point2f0(atan(0.5(x[2] - c[2]), 0.5(x[1] - c[1])), Q(x, p₁))
      end,
      markersize = 16px,
      color = :MEDIUMAQUAMARINE,
    )

    shadow_energy_plt = scatter!(
      energy_ax,
      lift(shadow_x, p₁) do x, p₁
        Point2f0(atan(0.5(x[2] - c[2]), 0.5(x[1] - c[1])), Q(x, p₁))
      end,
      markersize = 16px,
      marker = :diamond,
      color = :lightgrey,
    )
    shadow_energy_plt.visible[] = shadow_visible[]
    energy_ax.ytickformat = xs -> [@sprintf("%.3f", x) for x in xs]
    energy_ax.xticks = -π:π/2:π
    energy_ax.xtickformat = xs -> ["-π", "-π/2", "0", "π/2", "π"]

    autolimits!(energy_ax)
    on(energy_landscape) do _
      autolimits!(energy_ax)
    end

    catastrophe_btn = LButton(
      scene,
      label = "Show catastrophe curve",
      buttoncolor = :transparent,
      strokecolor = :black,
      tellwidth = false,
      width = 260
    )
    buttongrid[1, 1] = catastrophe_btn
    on(catastrophe_btn.clicks) do a
      active = isodd(a)
      catastrophe_curve.visible[] = active
      if active
        catastrophe_btn.label = "Hide catastrophe curve"
        shadow_visible[] = shadow_visible[]
      else
        catastrophe_btn.label = "Show catastrophe curve"
      end

    end
    catastrophe_curve.visible[] = false

    on(shadow_visible) do visible
      if isodd(catastrophe_btn.clicks[])
        shadow_energy_plt.visible[] = visible
        shadow_plt.visible[] = visible
      end
    end


    display(scene)
  end
end

end # module
