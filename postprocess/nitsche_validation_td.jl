# nitsche_validation_td.jl
#
# Validates enforcement of the Nitsche kinematic constraint for Case 2's
# time-domain formulation. Runs a short single-frequency simulation
# (n_periods_val wave periods, αp = 0 for the cleanest FSI check, no
# porous damping) and reports the constraint residual
#
#   J = u·n - ζ = 0   on Γstr,   where ζ = ∂t(η)
#
# normalised as ‖J‖/(ω η₀), over time. The mean residual over the
# analysis window [val_win_start T, val_win_end T] is the number quoted
# in the methodology chapter 4 as evidence that the constraint is
# well-enforced.
#
# Runs independently — does not require case2_sweep_L4.jl/L18.jl to
# have been run first.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "Thesis_FPVBarrier.jl"))
using .Thesis_FPVBarrier
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Plots
using JLD2
using Printf
using LaTeXStrings
import Statistics

# ── Settings ──────────────────────────────────────────────────────────────────
order          = 2
αp_val         = 0.0    # αp=0: no porous damping, cleanest FSI check
n_periods_val  = 15      # total periods to run
dt_per_T       = 40     # timesteps per period
mesh_file      = "meshes/case2_length_sweep/case2_L_18.msh"
val_win_start  = 12    # analysis window start [T]
val_win_end    = 15    # analysis window end   [T]
plot_dir       = "results/plots/validation/nitsche"
mkpath(plot_dir)

println("\n── Nitsche validation (time-domain) ────────────────────────────")
println("   mesh        : ", mesh_file)
println("   αp          : ", αp_val, " s⁻¹")
println("   n_periods   : ", n_periods_val)
println("   dt_per_T    : ", dt_per_T)
println("   window      : [$(val_win_start)T, $(val_win_end)T]")
println("────────────────────────────────────────────────────────────────")

# ── Build model ───────────────────────────────────────────────────────────────
model2  = Thesis_FPVBarrier.build_model(mesh_file)
params2 = Thesis_FPVBarrier.build_case2_transient_params(model2, αp_val)
reg2    = Thesis_FPVBarrier.build_case2_transient_regions(model2)
meas2   = Thesis_FPVBarrier.build_case2_transient_measures(reg2; order=order)
sp2     = Thesis_FPVBarrier.build_case2_transient_spaces(reg2, params2; order=order)
op2     = Thesis_FPVBarrier.build_case2_transient_operator(sp2, reg2, meas2, order, params2)

(; T, η₀, ω, g, H, k, Lb, γ₀, ρw) = params2

dt   = T / dt_per_T
tend = n_periods_val * T

# Mean element size on Γstr — for γN reporting
xΓstr = get_cell_coordinates(reg2.Γstr)
h_elm = sum(norm(xs[2] - xs[1]) for xs in xΓstr) / length(xΓstr)

println("  γ₀                = ", γ₀)
println("  γN = γ₀/h         = ", round(γ₀ / h_elm, digits=2))
println("  ρw·ω·h (inertial) = ", round(ρw * ω * h_elm, digits=2))
println("  kLb               = ", round(k * Lb, digits=3))
println()

# ── Initial conditions ────────────────────────────────────────────────────────
p0_fun = interpolate_everywhere(x -> g*(H - x[2]), sp2.U_p)
x₀     = interpolate_everywhere(
            [VectorValue(0.0, 0.0), p0_fun, 0.0, 0.0, 0.0], sp2.X(0.0));

ls         = LUSolver()
ode_solver = ThetaMethod(ls, dt, 0.5)

# ── Storage ───────────────────────────────────────────────────────────────────
J_nitsche_t = Float64[]
ts_val      = Float64[]

# ── Time loop ─────────────────────────────────────────────────────────────────
println("Starting validation loop...")
xht = solve(ode_solver, op2, x₀, 0.0, tend)

@time for ((uh, ph, κh, ηh, ζh), t) in xht

    # Kinematic residual: J = u·n - ζ  (ζ = ∂t(η) is the auxiliary DOF)
    J    = (uh ⋅ reg2.nstr) - ζh
    J_sq = sum( ∫( J * J ) * meas2.dΓstr )

    push!(J_nitsche_t, sqrt(abs(J_sq)) / (ω * η₀))
    push!(ts_val, t)

    # Progress
    if mod(length(ts_val), dt_per_T) == 0
        n_p = div(length(ts_val), dt_per_T)
        @printf("  t = %2d T  |  ‖J‖/(ω η₀) = %.3e\n", n_p, J_nitsche_t[end])
    end
end

println("\nValidation loop complete.")

# ── Window statistics ─────────────────────────────────────────────────────────
t_norm   = ts_val ./ T
mask     = (t_norm .>= val_win_start) .& (t_norm .<= val_win_end)
J_window = J_nitsche_t[mask]

println("\n── Results ─────────────────────────────────────────────────────")
if !isempty(J_window)
    J_mean = Statistics.mean(J_window)
    J_max  = Statistics.maximum(J_window)
    @printf("  Analysis window [%.1fT, %.1fT]:\n", val_win_start, val_win_end)
    @printf("    mean ‖J‖/(ω η₀) = %.4e\n", J_mean)
    @printf("    max  ‖J‖/(ω η₀) = %.4e\n", J_max)
    if J_mean < 0.01
        println("  ✓ Kinematic constraint well-satisfied (mean < 0.01)")
    elseif J_mean < 0.05
        println("  ~ Acceptable (mean < 0.05) — consider increasing γ₀")
    else
        println("  ✗ Poor enforcement — increase γ₀")
    end
end
println("────────────────────────────────────────────────────────────────")

# ── Plot ──────────────────────────────────────────────────────────────────────
p_J = plot(t_norm, J_nitsche_t,
    xlabel        = L"t \, / \, T",
    ylabel        = "Nitsche residual \$||J|| / (\\omega \\, \\eta_0)\$ ",
    # title         = "Nitsche kinematic residual",
    label         =  "γ₀ = $(γ₀)",
    lw            = 1.5,
    color         = :darkgreen,
    yscale        = :log10,
    yminorgrid    = true,
    framestyle    = :box,
    bottom_margin = 2Plots.mm,
    left_margin   = 2Plots.mm,
    right_margin  = 2Plots.mm,
    legend        =:topleft)

hline!(p_J, [0.01], linestyle=:dash, color=:gray,   label="target < 0.01")
hline!(p_J, [0.05], linestyle=:dot,  color=:orange, label="acceptable < 0.05")
# vline!(p_J, [val_win_start], linestyle=:dot, color=:black, lw=1.0, label="win start")
# vline!(p_J, [val_win_end],   linestyle=:dot, color=:gray,  lw=1.0, label="win end")

savefig(p_J, joinpath(plot_dir, "case2_nitsche_J_td.png"))
display(p_J)
println("Plot saved: case2_nitsche_J_td.png")

# ── Save ──────────────────────────────────────────────────────────────────────
jld2_path = joinpath("results/jld2/nitsche_validation", "case2_nitsche_validation_td.jld2")
@save jld2_path J_nitsche_t ts_val T ω η₀ val_win_start val_win_end
println("JLD2 saved: ", jld2_path)

println("\n── Nitsche validation complete ──────────────────────────────────")