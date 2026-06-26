# case2_khab_validation.jl
#
# Validates the Case 2 time-domain formulation against the
# Khabakhpasheva et al. benchmark: builds a spatially-varying EI
# profile and structural parameters matching the published reference
# geometry, runs the time integration to t = 80T, and plots snapshots
# of the beam deflection envelope against both the Khabakhpasheva and
# Riyansyah reference data.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include("../src/Thesis_FPVBarrier.jl")
using .Thesis_FPVBarrier
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.CellData: mean
using JLD2
using WriteVTK

# -- Khabakhpasheva params ------------------------------------------------
function build_case2_khabakhpasheva_params(model, αₚ = 0.0)
    # Read geometry from mesh
    p = Thesis_FPVBarrier.build_case2_transient_params(model, αₚ)

    (; Lb, xb₀, H) = p

    # Wave parameters
    η₀      = 0.01
    α_ratio = 0.249        # ratio λ/Lb defining the incident wavelength relative to beam length
    λ       = α_ratio * Lb
    k       = 2π / λ
    g       = 9.81
    ω       = sqrt(g * k * tanh(k * H))
    T       = 2π / ω

    # Structural parameters
    ρw  = 1025.0
    ρb  = 900.0
    hb  = 0.2
    d₀  = 8.1561e-3
    m   = ρw * d₀

    # Spatially varying EI
    EI₁    = 47100.0       # first 20% of beam
    EI₂    = 471.0         # remaining 80%
    EI_fun = x -> x[1] < (xb₀ + 0.2*Lb) ? EI₁ : EI₂
    EI     = EI₂           # scalar fallback (not used in khab operator)

    return merge(p, (
        η₀=η₀, λ=λ, k=k, ω=ω, T=T, kp=k, ωp=ω,
        ρw=ρw, ρb=ρb, hb=hb, d₀=d₀, m=m,
        EI=EI, EI_fun=EI_fun,
        EI₁=EI₁, EI₂=EI₂,
    ))
end

# -- Warmup ---------------------------------------------------------------------------
println("Warming up...")
let
    order    = 2
    model_wm = Thesis_FPVBarrier.build_model("meshes/validation/case2_validation_khab.msh")
    params_wm = build_case2_khabakhpasheva_params(model_wm, 0.0)
    reg_wm   = Thesis_FPVBarrier.build_case2_transient_regions(model_wm)
    meas_wm  = Thesis_FPVBarrier.build_case2_transient_measures(reg_wm; order=order)
    sp_wm    = Thesis_FPVBarrier.build_case2_transient_spaces(reg_wm, params_wm; order=order)
    op_wm    = Thesis_FPVBarrier.build_case2_khab_transient_operator(sp_wm, reg_wm, meas_wm, order, params_wm)
    p0_wm    = interpolate_everywhere(x -> params_wm.g*(params_wm.H - x[2]), sp_wm.U_p)
    x_wm     = interpolate_everywhere([VectorValue(0.0,0.0), p0_wm, 0.0, 0.0, 0.0], sp_wm.X(0.0))
    ls_wm    = LUSolver()
    xht_wm   = solve(ThetaMethod(ls_wm, params_wm.T/2, 0.5), op_wm, x_wm, 0.0, params_wm.T/2)
    for _ in xht_wm end
end
println("Warmup complete.")

# -- Main validation run ------------------------------------------------------------
function run_case2_khab_validation()
    order    = 2
    dt_per_T = 40
    tend_T   = 80

    mesh_path = "meshes/validation/case2_validation_khab.msh"
    jld2_path = "results/jld2/case2_validation/case2_khab_validation.jld2"
    vtk_dir   = "results/vtk/validation/case2"
    mkpath(dirname(jld2_path))
    mkpath(vtk_dir)

    model  = Thesis_FPVBarrier.build_model(mesh_path)
    params = build_case2_khabakhpasheva_params(model, 0.0)
    reg    = Thesis_FPVBarrier.build_case2_transient_regions(model)
    meas   = Thesis_FPVBarrier.build_case2_transient_measures(reg; order=order)
    sp     = Thesis_FPVBarrier.build_case2_transient_spaces(reg, params; order=order)
    op     = Thesis_FPVBarrier.build_case2_khab_transient_operator(sp, reg, meas, order, params)

    (; T, η₀, g, H, k, ω, LΩ, xb₀, xb₁, λ, m, ρw) = params

    dt   = T / dt_per_T
    tend = tend_T * T

    println("\n========================================")
    println("Case 2 Khabakhpasheva validation")
    println("  λ  = $(round(λ,  digits=4)) m")
    println("  k  = $(round(k,  digits=4)) rad/m")
    println("  ω  = $(round(ω,  digits=4)) rad/s")
    println("  T  = $(round(T,  digits=4)) s")
    println("  H  = $(round(H,  digits=4)) m")
    println("  LΩ = $(round(LΩ, digits=1)) m")
    println("  dt = T/$(dt_per_T) = $(round(dt, digits=5)) s")
    println("  tend = $(tend_T)T = $(round(tend, digits=3)) s")
    println("========================================\n")

    # -- Spatial DOF coordinates ------------------------------------------------------
    x_κ_all    = get_free_dof_values(interpolate_everywhere(x -> x[1], sp.V_κ))
    sort_idx_κ = sortperm(x_κ_all)
    xs_κ       = x_κ_all[sort_idx_κ]

    x_η_all    = get_free_dof_values(interpolate_everywhere(x -> x[1], sp.V_η))
    sort_idx_η = sortperm(x_η_all)
    xs_η       = x_η_all[sort_idx_η]

    # -- Storage ----------------------------------------------------------------------
    κ_profiles = Vector{Vector{Float64}}()
    η_profiles = Vector{Vector{Float64}}()
    ts         = Float64[]

    # -- Initial conditions -----------------------------------------------------------
    p0_fun = interpolate_everywhere(x -> g*(H - x[2]), sp.U_p)
    x₀     = interpolate_everywhere([VectorValue(0.0,0.0), p0_fun, 0.0, 0.0, 0.0], sp.X(0.0))

    ls         = LUSolver()
    ode_solver = ThetaMethod(ls, dt, 0.5)
    xht        = solve(ode_solver, op, x₀, 0.0, tend)

    step             = 0
    steps_per_period = round(Int, T/dt)

    # -- VTK collections ------------------------------------------------------------
    pvd_Ω    = paraview_collection(joinpath(vtk_dir, "case2_khab_omega"), append=false)
    pvd_Γfs  = paraview_collection(joinpath(vtk_dir, "case2_khab_fs"),    append=false)
    pvd_Γstr = paraview_collection(joinpath(vtk_dir, "case2_khab_beam"),  append=false)


    println("Starting time integration...")
    @time for ((uh, ph, κh, ηh, ζh), t) in xht
        step += 1
        push!(ts, t)

        κ_dofs = get_free_dof_values(κh)
        η_dofs = get_free_dof_values(ηh)

        push!(κ_profiles, copy(κ_dofs[sort_idx_κ] ./ η₀))
        push!(η_profiles, copy(η_dofs[sort_idx_η] ./ η₀))

        if mod(step, steps_per_period) == 0
            n = div(step, steps_per_period)
            println("  t = $(n)T  |  max|κ|/η₀ = $(round(maximum(abs.(κ_dofs))/η₀, digits=3))  |  max|η|/η₀ = $(round(maximum(abs.(η_dofs))/η₀, digits=3))")
        end

        # VTK output every period
        if mod(step, 4) == 0
            pvd_Ω[t]    = createvtk(reg.Ω,    joinpath(vtk_dir, "case2_khab_omega_$(step)"),
                            cellfields=["u"=>uh, "p"=>ph])
            pvd_Γfs[t]  = createvtk(reg.Γfs,  joinpath(vtk_dir, "case2_khab_fs_$(step)"),
                            cellfields=["kappa"=>κh])
            pvd_Γstr[t] = createvtk(reg.Γstr, joinpath(vtk_dir, "case2_khab_beam_$(step)"),
                            cellfields=["eta"=>ηh, "zeta"=>ζh])
        end
    end

    vtk_save(pvd_Ω)
    vtk_save(pvd_Γfs)
    vtk_save(pvd_Γstr)

    println("\nTime integration complete. Saving JLD2...")
    jldopen(jld2_path, "w") do f
        f["κ_profiles"] = κ_profiles
        f["η_profiles"] = η_profiles
        f["xs_κ"]       = xs_κ
        f["xs_η"]       = xs_η
        f["ts"]         = ts
        f["T"]          = T
        f["ω"]          = ω
        f["k"]          = k
        f["η₀"]         = η₀
        f["LΩ"]         = LΩ
        f["xb₀"]        = xb₀
        f["xb₁"]        = xb₁
        f["H"]          = H
        f["λ"]          = λ
        f["m"]          = m
        f["ρw"]          = ρw
    end
    println("Saved: ", jld2_path)
end

@time run_case2_khab_validation()

# -- Postprocess and plotting  ----------------------------------------------------
using JLD2, Plots, LaTeXStrings, CSV
using Plots: mm

plot_dir = "results/plots/validation/case2"
mkpath(plot_dir)

default(fontfamily="DejaVu Sans", guidefontsize=11, tickfontsize=9,
        legendfontsize=8, titlefontsize=14)

# -- Load JLD2 -------------------------------------------------------------------
f          = load("results/jld2/case2_validation/case2_khab_validation.jld2")
T          = f["T"]
ts         = f["ts"]
xs_η       = f["xs_η"]
η_profiles = f["η_profiles"]
LΩ         = f["LΩ"]
xb₀        = f["xb₀"]
xb₁        = f["xb₁"]
Lb         = xb₁ - xb₀

# -- Reference data ----------------------------------------------------------------
ref_khab = CSV.File("results/data/ref_data/Khabakhpasheva_without_joint.csv"; header=false)
ref_riy  = CSV.File("results/data/ref_data/Riyansyah_without_joint.csv";      header=false)

# -- Snapshot times ----------------------------------------------------------------
# Four snapshots within roughly one period, spaced 0.2T apart, centred on
# the peak beam deflection observed during the run.
t_max = 47.0
snap_Ts = [t_max, t_max + 0.2, t_max + 0.4, t_max + 0.6, t_max + 0.8]
colors_snap = cgrad(:blues, length(snap_Ts), categorical=true)

# -- Plot --------------------------------------------------------------------------
plt_η = plot(
    xlabel     = L"x / L_b",
    ylabel     = L"|\eta| / \eta_0",
    # title      = "Beam deflection envelope, (\$\\alpha_p = 0\$)",
    framestyle = :box, grid = true,
    size       = (1000, 350),
    top_margin    = 2Plots.mm, bottom_margin = 5Plots.mm,
    left_margin   = 5Plots.mm, right_margin  = 5Plots.mm,
    legend        = :outerright
)

for (i, t_T) in enumerate(snap_Ts)
    i_snap  = argmin(abs.(ts .- t_T * T))
    η_snap  = abs.(η_profiles[i_snap])
    x_norm  = (xs_η .- xb₀) ./ Lb
    plot!(plt_η, x_norm, η_snap,
        label = "t = $(round(t_T, digits=2))T",
        lw    = 1.5,
        color = colors_snap[i])
end

# Reference data
plot!(plt_η,
    collect(ref_khab.Column1), collect(ref_khab.Column2),
    seriestype = :scatter, marker = :circle, ms = 4,
    color = :black, label = "Khabakhpasheva et al.")

plot!(plt_η,
    collect(ref_riy.Column1), collect(ref_riy.Column2),
    lw = 2, ls = :dash, color = :red, label = "Riyansyah et al.")

vline!(plt_η, [0.0, 1.0], ls = :dot, color = :blue, lw = 1.0, label = false)

savefig(plt_η, joinpath(plot_dir, "case2_khab_eta_envelope.png"))
display(plt_η)
println("Saved: case2_khab_eta_envelope.png")