# case3_postprocess.jl
#
# Post-processing for Case 3: floating rigid pipe, frequency domain.
# Requires case3.jl to have been run first (provides params3, reg3,
# meas3, sp3, κh3, ξ).
#
# Runs the full sequence of analyses behind every Case 3 figure and
# table in the thesis:
#
#   1-4. Single-frequency spatial plots and rigid-body DOFs, at the
#        base case already solved by case3.jl
#   5.   Base frequency sweep (D = 0.45 m, unmoored)
#   6.   Diameter sweep, unmoored, D ∈ {0.25, ..., 0.50} m
#   7.   Coarse mooring stiffness sweep, C = 0, at ωp (50 log-spaced
#        points) — identifies k/Khs11 ≈ 1.151 as the approximate
#        optimum, common to all diameters at this resolution
#   8.   Refined mooring stiffness sweep (50 linearly-spaced points
#        around 1.151) — resolves the per-diameter variation in the
#        true optimum that the coarse sweep could not detect
#   9.   Moored diameter sweep, using the refined optimal stiffness
#        per diameter
#  10.   Willemspolder as-built mooring stiffness, derived from the
#        mooring rope geometry and compared against the unmoored and
#        optimally-moored cases
#
# All sweep results are saved to a single JLD2 file
# (results/jld2/case3/case3_results.jld2) so they can be reloaded for
# further analysis or replotting without rerunning the sweeps — see
# load_case3_results() near the end of the file.

using Plots
using Printf
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces
using LaTeXStrings
using JLD2

plot_dir = "results/plots/case3"
mkpath(plot_dir)
jld2_dir_case3 = "results/jld2/case3"
mkpath(jld2_dir_case3)

default(fontfamily="DejaVu Sans", guidefontsize=11, tickfontsize=9, legendfontsize=8, titlefontsize=14)

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: extract sorted (x, κ) from κh using reg.Γfs directly
# ─────────────────────────────────────────────────────────────────────────────
function extract_kappa_case3(κh, reg)
    xΓfs       = get_cell_coordinates(reg.Γfs)
    κ_dof_vals = get_cell_dof_values(κh)
    x_vals = Float64[]
    κ_vals = ComplexF64[]
    for i in 1:num_cells(reg.Γfs)
        cell_coords = xΓfs[i]
        xm   = sum(cell_coords[j][1] for j in 1:length(cell_coords)) / length(cell_coords)
        dofs = κ_dof_vals[i]
        push!(x_vals, xm)
        push!(κ_vals, sum(dofs) / length(dofs))
    end
    idx = sortperm(x_vals)
    return x_vals[idx], κ_vals[idx]
end

# ─────────────────────────────────────────────────────────────────────────────
# 1. COMBINED SPATIAL PLOT: |κ|, |κ_r|, |κ_in| vs x/LΩ
#    κ split at pipe waterline gap, pipe location shown as shaded region
# ─────────────────────────────────────────────────────────────────────────────
function plot_kappa_eta_case3(κh, reg, params, ξ)
    (; η₀, LΩ, k, xdᵢₙ, xdₒᵤₜ, xc, xpL, xpR) = params

    x_fs, κ_vals = extract_kappa_case3(κh, reg)
    κ_in_vals    = η₀ .* exp.(im .* k .* x_fs)
    κ_r_vals     = κ_vals .- κ_in_vals

    κ_abs    = abs.(κ_vals)    ./ η₀
    κ_r_abs  = abs.(κ_r_vals)  ./ η₀
    κ_in_abs = abs.(κ_in_vals) ./ η₀

    # Pipe heave motion, reported as a single marker at the pipe centre,
    # since the rigid body has no spatial deflection profile to plot —
    # only the overall heave/pitch response
    η_p_heave = abs(ξ[1]) / η₀

    mask_left  = x_fs .< xpL
    mask_right = x_fs .> xpR

    # Plot 1: full domain
    p1 = plot(
        x_fs[mask_left] ./ LΩ, κ_abs[mask_left],
        label="Free Surface", lw=2, color=:black,
        xlabel=L"x / L_{\Omega}", ylabel=L"Normalised amplitude $|\cdot| / \eta_0$",
        # title="Surface elevation amplitudes — Case 3",
        legend=:topright
    )
    plot!(p1, x_fs[mask_right] ./ LΩ, κ_abs[mask_right],   label="", lw=2, color=:black)
    plot!(p1, x_fs[mask_left]  ./ LΩ, κ_r_abs[mask_left],  label="Reflected wave", lw=1.5, color=:royalblue)
    plot!(p1, x_fs[mask_left]  ./ LΩ, κ_in_abs[mask_left], label="Incident wave", lw=1.0, color=:gray, ls=:dot)

    vline!(p1, [xdᵢₙ  / LΩ], color=:green, ls=:dash, lw=0.8, label="Damping edges")
    vline!(p1, [xdₒᵤₜ / LΩ], color=:green, ls=:dash, lw=0.8, label="")
    vline!(p1, [xc    / LΩ], color=:black,  ls=:dot,  lw=0.8, label="Pipe centre")
    vspan!(p1, [xpL/LΩ, xpR/LΩ], alpha=0.15, color=:red, label="Pipe waterline gap")

    display(p1)
    savefig(p1, joinpath(plot_dir, "case3_kappa_eta_spatial.png"))

    # Plot 2: Without damping zone
    small_mask_left  = (x_fs .>= xdᵢₙ) .& (x_fs .< xpL)
    small_mask_right = (x_fs .> xpR)   .& (x_fs .<= xdₒᵤₜ)

    p2 = plot(
        x_fs[small_mask_left] ./ LΩ, κ_abs[small_mask_left],
        label="Free Surface", lw=2, color=:black,
        xlabel=L"x / L_{\Omega}", ylabel=L"Normalised amplitude $|\cdot| / \eta_0$",
        # title="Surface elevation amplitudes — Case 3"
    )
    plot!(p2, x_fs[small_mask_right] ./ LΩ, κ_abs[small_mask_right],   label="", lw=2, color=:black)
    plot!(p2, x_fs[small_mask_left]  ./ LΩ, κ_r_abs[small_mask_left],  label="Reflected wave", lw=1.5, color=:royalblue)
    plot!(p2, x_fs[small_mask_left]  ./ LΩ, κ_in_abs[small_mask_left], label="Incident wave", lw=1.0, color=:gray, ls=:dot)
    vline!(p2, [xc / LΩ], color=:black, ls=:dot, lw=0.8, label="Pipe centre")
    vspan!(p2, [xpL/LΩ, xpR/LΩ], alpha=0.15, color=:red, label="Pipe waterline gap")
    scatter!(p2, [xc / LΩ], [η_p_heave],
            label=L"Pipe heave $|\xi_1|/\eta_0$", color=:red, ms=6, markerstrokewidth=0.5, markerstrokecolor=:black)

    display(p2)
    savefig(p2, joinpath(plot_dir, "case3_kappa_eta_spatial_inner.png"))
    return p1, p2
end

# ─────────────────────────────────────────────────────────────────────────────
# 2.  |κ|/η₀ along the free surface
# ─────────────────────────────────────────────────────────────────────────────
function plot_kappa_free_surface_case3(κh, reg, params)
    (; η₀, LΩ, xdᵢₙ, xdₒᵤₜ, xc) = params

    x_s, κ_vals = extract_kappa_case3(κh, reg)
    κ_norm = abs.(κ_vals) ./ η₀

    p = plot(
        x_s ./ LΩ, κ_norm,
        xlabel = L"x / L_{\Omega}", ylabel = L"$|\kappa| / \eta_0$",
        title  = L"$|\kappa| / \eta_0$ along free surface — Case 3",
        label  = L"$|\kappa| / \eta_0$", lw = 1.5
    )
    vline!(p, [xdᵢₙ  / LΩ], ls=:dash, color=:red,   label="Damping zone edge (in)")
    vline!(p, [xdₒᵤₜ / LΩ], ls=:dash, color=:green, label="Damping zone edge (out)")
    vline!(p, [xc    / LΩ], ls=:dot,  color=:black,  label="Pipe centre")

    display(p)
    savefig(p, joinpath(plot_dir, "case3_kappa_abs.png"))
    return x_s, κ_norm
end

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Re and Im of κ
# ─────────────────────────────────────────────────────────────────────────────
function plot_kappa_re_im_case3(κh, reg, params)
    (; η₀, LΩ, xdᵢₙ, xdₒᵤₜ, xpL, xpR, k) = params

    x_s, κ_s = extract_kappa_case3(κh, reg)
    κ_in = η₀ .* exp.(im .* k .* x_s)

    p = plot(
        x_s ./ LΩ, real.(κ_s) ./ η₀,
        xlabel = L"$x / L_{\Omega}$", ylabel = L"$\kappa / \eta_0$",
        title  = L"$Re$ and $Im$ of $\kappa$ — Case 3",
        label  = L"$Re(\kappa)/\eta_0$", lw = 1.5, color = :royalblue
    )
    plot!(p, x_s ./ LΩ, imag.(κ_s) ./ η₀,
        label = L"$Im(\kappa)/\eta_0$", lw = 1.5, color = :coral, ls = :dash)
    plot!(p, x_s ./ LΩ, real.(κ_in) ./ η₀,
        label = L"$Re(\kappa_{in})/\eta_0$", lw = 1.0, color = :gray, ls = :dot)
    vline!(p, [xdᵢₙ/LΩ, xdₒᵤₜ/LΩ], color=:black, ls=:dash, lw=0.8, label="Damping edges")

    display(p)
    savefig(p, joinpath(plot_dir, "case3_kappa_re_im.png"))
    return x_s, κ_s
end

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Rigid body DOFs — print to console
# ─────────────────────────────────────────────────────────────────────────────
function print_rigid_body_dofs(ξ, params)
    (; η₀, k) = params
    println("─── Rigid body DOFs ──────────────────────────")
    println("  ξ_heave = ", ξ[1])
    println("  ξ_pitch = ", ξ[2])
    println("  RAO heave  |ξ₁|/η₀        = ", round(abs(ξ[1])/η₀,        digits=4))
    println("  RAO pitch  |ξ₂|/(k⋅η₀)    = ", round(abs(ξ[2])/(k*η₀),    digits=4))
    println("  Phase heave [deg]          = ", round(rad2deg(angle(ξ[1])), digits=2))
    println("  Phase pitch [deg]          = ", round(rad2deg(angle(ξ[2])), digits=2))
end

# ─── Run single-frequency plots ───────────────────────────────────────────────
plot_kappa_eta_case3(κh3, reg3, params3, ξ)
plot_kappa_free_surface_case3(κh3, reg3, params3)
plot_kappa_re_im_case3(κh3, reg3, params3)
print_rigid_body_dofs(ξ, params3)

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: Mansard-Funke 3-probe decomposition
# ─────────────────────────────────────────────────────────────────────────────
function decompose_wave(x_cells, κ_cells, k, xp1, xp2, xp3)
    i1 = argmin(abs.(x_cells .- xp1))
    i2 = argmin(abs.(x_cells .- xp2))
    i3 = argmin(abs.(x_cells .- xp3))
    xs = [x_cells[i1], x_cells[i2], x_cells[i3]]
    κs = [κ_cells[i1], κ_cells[i2], κ_cells[i3]]
    M  = [exp(im*k*xs[1])  exp(-im*k*xs[1])
          exp(im*k*xs[2])  exp(-im*k*xs[2])
          exp(im*k*xs[3])  exp(-im*k*xs[3])]
    AB = M \ κs
    return AB[1], AB[2]
end

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: extract κ cells from FEFunction
# ─────────────────────────────────────────────────────────────────────────────
function extract_kappa_cells(κh, reg)
    xΓfs_i  = get_cell_coordinates(reg.Γfs)
    κ_dof_i = get_cell_dof_values(κh)
    xv = Float64[]; κv = ComplexF64[]
    for j in 1:num_cells(reg.Γfs)
        coords = xΓfs_i[j]
        xm = sum(coords[l][1] for l in 1:length(coords)) / length(coords)
        push!(xv, xm)
        push!(κv, sum(κ_dof_i[j]) / length(κ_dof_i[j]))
    end
    idx = sortperm(xv)
    return xv[idx], κv[idx]
end

# ─────────────────────────────────────────────────────────────────────────────
# CORE: single frequency solve — returns K_R, K_T, RAOs, drift, PErr
# ─────────────────────────────────────────────────────────────────────────────
function solve_single_frequency(k, params_base, sp, reg, meas; params_override=NamedTuple())
    (; g, H, η₀, xc, xdᵢₙ, xdₒᵤₜ) = params_base

    ω = sqrt(g * k * tanh(k * H))
    λ = 2π / k
    T = 2π / ω

    params_i = merge(params_base, (k=k, ω=ω, λ=λ, T=T))
    params_i = merge(params_i, params_override)

    sys_i = Thesis_FPVBarrier.build_case3_operator(sp, reg, meas, 2, params_i)
    x_i   = sys_i.A \ sys_i.b

    nfe_i = num_free_dofs(sp.X)
    xfe_i = x_i[1:nfe_i]
    ξ_i   = x_i[nfe_i+1:end]

    nϕ_i = num_free_dofs(sp.U_Ω)
    nκ_i = num_free_dofs(sp.U_κ)
    κh_i = FEFunction(sp.U_κ, xfe_i[nϕ_i+1:nϕ_i+nκ_i])

    x_cells, κ_cells = extract_kappa_cells(κh_i, reg)

    # Upstream 3-probe
    x_up_end = xc - (xc - xdᵢₙ) * 0.05
    x_p1     = xdᵢₙ + (xc - xdᵢₙ) * 0.30
    Δx_up    = min(λ / 4, (x_up_end - x_p1) / 2.5)
    A_inc, B_ref = decompose_wave(x_cells, κ_cells, k,
                       x_p1, x_p1 + Δx_up, x_p1 + 2Δx_up)
    K_R = abs(B_ref) / abs(A_inc)

    # Downstream 3-probe
    x_dn_start = xc     + (xdₒᵤₜ - xc) * 0.10
    x_dn_end   = xdₒᵤₜ - (xdₒᵤₜ - xc) * 0.05
    Δx_dn      = min(λ / 4, (x_dn_end - x_dn_start) / 2.5)
    C_trans, _ = decompose_wave(x_cells, κ_cells, k,
                       x_dn_start, x_dn_start + Δx_dn, x_dn_start + 2Δx_dn)
    K_T = abs(C_trans) / abs(A_inc)

    F_drift_nd = 1.0 + K_R - K_T
    PErr       = (1.0 - K_R^2 - K_T^2) * 100
    rao_h      = abs(ξ_i[1]) / η₀
    rao_p      = abs(ξ_i[2]) / (k * η₀)

    return (ω=ω, K_R=K_R, K_T=K_T, F_drift=F_drift_nd, PErr=PErr,
        rao_h=rao_h, rao_p=rao_p,
        phase_h=rad2deg(angle(ξ_i[1])), phase_p=rad2deg(angle(ξ_i[2])),
        Khs=sys_i.Khs,
        x_p1=x_p1, x_p2=x_p1+Δx_up, x_p3=x_p1+2Δx_up,
        x_dn1=x_dn_start, x_dn2=x_dn_start+Δx_dn, x_dn3=x_dn_start+2Δx_dn,
        Δx_up=Δx_up, λ=λ)
end

# ─────────────────────────────────────────────────────────────────────────────
# CORE: full frequency sweep — returns named tuple of result arrays
# ─────────────────────────────────────────────────────────────────────────────
function run_frequency_sweep(k_range, params_base, sp, reg, meas;
                              params_override=NamedTuple(), label="")
    N = length(k_range)
    ω_vals      = Float64[]
    K_R_vals    = Float64[]
    K_T_vals    = Float64[]
    RAO_heave   = Float64[]
    RAO_pitch   = Float64[]
    phase_heave = Float64[]
    phase_pitch = Float64[]
    F_drift     = Float64[]
    PErr_vals   = Float64[]

    println("Frequency sweep: $label ($(N) frequencies)...")
    flush(stdout)

    for (i, k) in enumerate(k_range)
        r = solve_single_frequency(k, params_base, sp, reg, meas;
                                   params_override=params_override)
        push!(ω_vals,      r.ω)
        push!(K_R_vals,    r.K_R)
        push!(K_T_vals,    r.K_T)
        push!(RAO_heave,   r.rao_h)
        push!(RAO_pitch,   r.rao_p)
        push!(phase_heave, r.phase_h)
        push!(phase_pitch, r.phase_p)
        push!(F_drift,     r.F_drift)
        push!(PErr_vals,   r.PErr)

        if mod(i, 10) == 0 || i == 1
            @printf("  [%2d/%d]  k=%.3f  ω=%.3f  K_R=%.3f  K_T=%.3f  PErr=%.1f%%\n",
                i, N, k, r.ω, r.K_R, r.K_T, r.PErr)
            flush(stdout)
        end
    end

    println("  Done: $label")
    return (ω=ω_vals, K_R=K_R_vals, K_T=K_T_vals,
            RAO_heave=RAO_heave, RAO_pitch=RAO_pitch,
            phase_heave=phase_heave, phase_pitch=phase_pitch,
            F_drift=F_drift, PErr=PErr_vals)
end

# ─────────────────────────────────────────────────────────────────────────────
# 5. BASE FREQUENCY SWEEP (D = 0.45 m, no mooring)
# ─────────────────────────────────────────────────────────────────────────────
(; g, H, η₀, ρw, xc, xdᵢₙ, xdₒᵤₜ, D, ωp, kp) = params3

# Frequency sweep ranges
ω_low  = range(0.5*ωp, 2.0, length=30)
ω_mid  = range(2.0, ωp, length=16)
ω_high = range(ωp, 2.0*ωp, length=16)
ω_sweep = unique(vcat(collect(ω_low), collect(ω_mid), collect(ω_high)))

# Convert ω to k via dispersion relation
function find_k(ω, H, g=9.81; tol=1e-10, maxiter=100)
    k = ω^2 / g  # deep water initial guess
    for _ in 1:maxiter
        f  = g * k * tanh(k * H) - ω^2
        df = g * (tanh(k * H) + k * H * (1 - tanh(k * H)^2))
        k  = k - f / df
        abs(f) < tol && break
    end
    return k
end

k_sweep = find_k.(ω_sweep, H)

@time res_base = run_frequency_sweep(k_sweep, params3, sp3, reg3, meas3; label="base D=$(D)m")

# Single frequency spatial plot at lowest frequency
k_low    = k_sweep[1]
ω_low_val = ω_sweep[1]
λ_low    = 2π / k_low
params_low = merge(params3, (k=k_low, ω=ω_low_val, λ=λ_low, T=2π/ω_low_val))
sys_low  = Thesis_FPVBarrier.build_case3_operator(sp3, reg3, meas3, 2, params_low)
x_low    = sys_low.A \ sys_low.b
nfe_low  = num_free_dofs(sp3.X)
nϕ_low   = num_free_dofs(sp3.U_Ω)
nκ_low   = num_free_dofs(sp3.U_κ)
κh_low   = FEFunction(sp3.U_κ, x_low[1:nfe_low][nϕ_low+1:nϕ_low+nκ_low])

plot_kappa_free_surface_case3(κh_low, reg3, params_low)


# Panel plot
p_KR = plot(res_base.ω, res_base.K_R,
    xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel=L"K_{R}", title="Reflection coefficient",
    label=L"K_{R}", lw=2, color=:royalblue, ylims=(0, 1.1))
vline!(p_KR, [ωp], color=:gray, ls=:dash, lw=0.8, label=L"\omega_{p}")

p_KT = plot(res_base.ω, res_base.K_T,
    xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel=L"K_{T}", title="Transmission coefficient",
    label=L"K_{T}", lw=2, color=:coral, ylims=(0, 1.1))
vline!(p_KT, [ωp], color=:gray, ls=:dash, lw=0.8, label=L"\omega_{p}")

p_drift = plot(res_base.ω, res_base.F_drift,
    xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel=L"F_d / (IC_g)_\mathrm{in}", title="Drift force",
    label=L"1 + K_{R} - K_{T}", lw=2, color=:purple)
vline!(p_drift, [ωp], color=:gray, ls=:dash, lw=0.8, label=L"\omega_{p}")

p_err = plot(res_base.ω, res_base.PErr,
    xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel=L"PErr \ [\%]", title="Energy conservation error",
    label=L"1 - K_{R}^2 - K_{T}^2", lw=2, color=:gray)
hline!(p_err, [0], color=:black, ls=:dash, lw=0.8, label="")

p_panel = plot(p_KR, p_KT, p_drift, p_err,
    layout=(2,2), size=(900,700),
    plot_title="Case 3 — base frequency sweep (D=$(D)m)")
display(p_panel)
savefig(p_panel, joinpath(plot_dir, "case3_frequency_sweep_base.png"))

p_rao = plot(res_base.ω, res_base.RAO_heave,
    xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel="RAO", title="Response amplitude operators",
    label=L"Heave $|\xi_1| / \eta_0$", lw=2, color=:royalblue)
plot!(p_rao, res_base.ω, res_base.RAO_pitch,
    label=L"Pitch $|\xi_2| / (k \eta_0)$", lw=2, color=:coral, ls=:dash)
hline!(p_rao, [1.0], color=:gray, ls=:dot, lw=0.8, label="RAO = 1")
vline!(p_rao, [ωp], color=:gray, ls=:dash, lw=0.8, label=L"\omega_{p}")
display(p_rao)
savefig(p_rao, joinpath(plot_dir, "case3_RAOs_base.png"))



# ─────────────────────────────────────────────────────────────────────────────
# 6. DIAMETER SWEEP
# ─────────────────────────────────────────────────────────────────────────────
D_vals   = [0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
D_colors = [:royalblue, :cyan, :green, :orange, :red, :purple]
D_meshes = [
    "meshes/case3_diametersweep/case3_D_25.msh",
    "meshes/case3_diametersweep/case3_D_30.msh",
    "meshes/case3_diametersweep/case3_D_35.msh",
    "meshes/case3_diametersweep/case3_D_40.msh",
    "meshes/case3_diametersweep/case3_D_45.msh",
    "meshes/case3_diametersweep/case3_D_50.msh",
]

# ── Run unmoored sweep and store results ─────────────────────────────────────
unmoored_results = []

@time for (D_i, col, msh_path) in zip(D_vals, D_colors, D_meshes)
    model_i  = Thesis_FPVBarrier.build_model(msh_path)
    params_i = Thesis_FPVBarrier.build_case3_params(model_i)
    reg_i    = Thesis_FPVBarrier.build_case3_regions(model_i, params_i)
    meas_i   = Thesis_FPVBarrier.build_case3_measures(reg_i; order=2)
    sp_i     = Thesis_FPVBarrier.build_case3_spaces(reg_i, params_i; order=2)

    res_u = run_frequency_sweep(k_sweep, params_i, sp_i, reg_i, meas_i;
                                label="D=$(D_i)m unmoored")

    # Get Khs11 and natural frequency for this diameter
    r_ref_i  = solve_single_frequency(kp, params_i, sp_i, reg_i, meas_i)
    Khs11_i  = real(r_ref_i.Khs[1,1])
    i_nat_i  = argmax(res_u.RAO_heave)
    ωn_i     = res_u.ω[i_nat_i]
    m_eff_i  = Khs11_i / ωn_i^2
    k_moor_opt_i = max(0.0, ωp^2 * m_eff_i - Khs11_i)

    # Sanity check natural frequency
    m_i      = real(params_i.Mrb[1,1])
    R_i      = params_i.Rₒᵤₜ
    m_add_i  = 1000 * π * R_i^2
    ωn_dry_i = sqrt(Khs11_i / m_i)
    ωn_wet_i = sqrt(Khs11_i / (m_i + m_add_i))
    @printf("D = %.2f m  ωn_dry = %.3f  ωn_wet = %.3f  ωn_RAO = %.3f\n",
        D_i, ωn_dry_i, ωn_wet_i, ωn_i)

    println("D = $(D_i) m  →  ωn = $(round(ωn_i, digits=3))  Khs11 = $(round(Khs11_i, digits=1))  m_eff = $(round(m_eff_i, digits=1))  k_moor_opt = $(round(k_moor_opt_i, digits=1))")

    push!(unmoored_results, (D=D_i, col=col, res=res_u,
                             Khs11=Khs11_i, ωn=ωn_i, m_eff=m_eff_i,
                             k_moor_opt=k_moor_opt_i,
                             m_i=m_i, m_add_i=m_add_i, R_i=R_i, ωn_dry_i=ωn_dry_i, ωn_wet_i=ωn_wet_i,
                             params=params_i, sp=sp_i, reg=reg_i, meas=meas_i))
end


p_peak = plot(xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel="RAO heave", 
              title="Peak region — small D comparison",
              xlims=(1.2, 2.0))

for entry in unmoored_results
    entry.D > 0.35 && continue
    mask = (entry.res.ω .>= 1.2) .& (entry.res.ω .<= 2.0)
    plot!(p_peak, entry.res.ω[mask], entry.res.RAO_heave[mask],
          label=L"D = %$(entry.D) \ \mathrm{m}", lw=2, color=entry.col,
          marker=:circle, markersize=3)
    vline!(p_peak, [entry.ωn], color=entry.col, ls=:dash, lw=0.8, label=L"ωn \, D=%$(entry.D)")
end

display(p_peak)

# ── Plot unmoored results ─────────────────────────────────────────────────────
p_KR_D   = plot(xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel=L"K_{R}", title="Reflection", ylims=(0,1.1))
p_KT_D   = plot(xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel=L"K_{T}", title="Transmission", ylims=(0,1.1), legend=:left)
p_rao_D  = plot(xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel="RAO heave", title="Heave RAO")

vline!(p_KR_D,  [ωp], color=:gray, ls=:dash, lw=0.8, label=L"\omega_{p}")
vline!(p_KT_D,  [ωp], color=:gray, ls=:dash, lw=0.8, label=L"\omega_{p}")
vline!(p_rao_D, [ωp], color=:gray, ls=:dash, lw=0.8, label=L"\omega_{p}")

for entry in unmoored_results
    lbl = L"D = %$(entry.D) \ \mathrm{m}" * (entry.D == 0.45 ? " ★" : "")
    lw  = entry.D == 0.45 ? 2.5 : 1.5

    plot!(p_KR_D,  entry.res.ω, entry.res.K_R,       label=lbl, lw=lw, color=entry.col)
    plot!(p_KT_D,  entry.res.ω, entry.res.K_T,       label=lbl, lw=lw, color=entry.col)
    plot!(p_rao_D, entry.res.ω, entry.res.RAO_heave, label=lbl, lw=lw, color=entry.col)
end

p_D_panel = plot(p_KR_D, p_KT_D, p_rao_D,
    layout=(1,3), size=(1200,400),
    top_margin=7Plots.mm,
    bottom_margin=7Plots.mm,
    left_margin=7Plots.mm,
    right_margin=7Plots.mm)
display(p_D_panel)
savefig(p_D_panel, joinpath(plot_dir, "case3_D_sweep_unmoored.png"))

# # ─────────────────────────────────────────────────────────────────────────────
# # 7. MOORING STIFFNESS SWEEP (per diameter, C = 0, at ωp)
# # ─────────────────────────────────────────────────────────────────────────────

k_heave_factors = exp10.(range(-3, 2, length=50))

kmoor_sweep_results = []

println("\nMooring stiffness sweep (per diameter)...")
@time for entry in unmoored_results
    r_ref   = solve_single_frequency(kp, entry.params, entry.sp, entry.reg, entry.meas)
    Khs11_i = real(r_ref.Khs[1,1])
    k_heave_vals = k_heave_factors .* Khs11_i

    KR_kmoor = Float64[]
    KT_kmoor = Float64[]

    for k_heave in k_heave_vals
        r = solve_single_frequency(kp, entry.params, entry.sp, entry.reg, entry.meas;
                params_override=(
                    Krb = ComplexF64[k_heave 0.0; 0.0 0.0],
                    Crb = ComplexF64[0.0     0.0; 0.0 0.0]
                ))
        push!(KR_kmoor, r.K_R)
        push!(KT_kmoor, r.K_T)
    end

    i_opt      = argmin(KT_kmoor)
    k_moor_opt = k_heave_vals[i_opt]
    @printf("D = %.2f m  k_opt/Khs11 = %.3f  K_T_min = %.3f\n",
        entry.D, k_heave_factors[i_opt], KT_kmoor[i_opt])

    push!(kmoor_sweep_results, (
        D          = entry.D,
        col        = entry.col,
        Khs11      = Khs11_i,
        k_moor_opt = k_moor_opt,
        factors    = k_heave_factors,
        KR_kmoor   = KR_kmoor,
        KT_kmoor   = KT_kmoor
    ))
end

# ── Plot mooring stiffness sweep per diameter ─────────────────────────────────
p_kmoor = plot(xlabel="k_heave / Khs11", ylabel="K_R, K_T",
    title="Mooring stiffness sweep, ω = $(round(ωp, digits=3)) rad/s",
    xscale=:log10, ylims=(0,1.1), legend=:left)

for entry in kmoor_sweep_results
    lbl = "D = $(entry.D) m"
    lw  = entry.D == 0.45 ? 2.5 : 1.5
    plot!(p_kmoor, entry.factors, entry.KR_kmoor,
        label=lbl*" K_R", lw=lw, color=entry.col)
    plot!(p_kmoor, entry.factors, entry.KT_kmoor,
        label=lbl*" K_T", lw=lw, color=entry.col, ls=:dash)
end

vline!(p_kmoor, [1.0], color=:gray, ls=:dot, lw=0.8, label="k=Khs11")
display(p_kmoor)
savefig(p_kmoor, joinpath(plot_dir, "case3_mooring_stiffness_sweep.png"))


# ─────────────────────────────────────────────────────────────────────────────
# 8. REFINED MOORING STIFFNESS SWEEP around k_opt/Khs11 = 1.151
# ─────────────────────────────────────────────────────────────────────────────
k_refine_factors = range(0.9, 1.3, length=50)  # fine linear sweep around 1.151

kmoor_refined_results = []

println("\nRefined mooring stiffness sweep...")
@time for entry in unmoored_results
    r_ref   = solve_single_frequency(kp, entry.params, entry.sp, entry.reg, entry.meas)
    Khs11_i = real(r_ref.Khs[1,1])
    k_heave_vals = k_refine_factors .* Khs11_i

    KT_kmoor = Float64[]
    KR_kmoor = Float64[]

    for k_heave in k_heave_vals
        r = solve_single_frequency(kp, entry.params, entry.sp, entry.reg, entry.meas;
                params_override=(
                    Krb = ComplexF64[k_heave 0.0; 0.0 0.0],
                    Crb = ComplexF64[0.0     0.0; 0.0 0.0]
                ))
        push!(KT_kmoor, r.K_T)
        push!(KR_kmoor, r.K_R)
    end

    i_opt      = argmin(KT_kmoor)
    k_moor_opt = k_heave_vals[i_opt]
    @printf("D = %.2f m  k_opt/Khs11 = %.4f  K_T_min = %.4f  K_R = %.4f  PErr = %.3f%%\n",
        entry.D, k_refine_factors[i_opt], KT_kmoor[i_opt], KR_kmoor[i_opt],
        (1 - KR_kmoor[i_opt]^2 - KT_kmoor[i_opt]^2) * 100)

    push!(kmoor_refined_results, (
        D          = entry.D,
        col        = entry.col,
        Khs11      = Khs11_i,
        k_moor_opt = k_moor_opt,
        factors    = k_refine_factors,
        KR_kmoor   = KR_kmoor,
        KT_kmoor   = KT_kmoor
    ))
end

# ── Plot refined sweep ────────────────────────────────────────────────────────
p_kmoor_ref = plot(xlabel="k_heave / Khs11", ylabel="K_T",
    # title="Refined mooring stiffness sweep, ω = $(round(ωp, digits=3)) rad/s",
    ylims=(0, 1.1), legend=:topright)

for entry in kmoor_refined_results
    lbl = "D = $(entry.D) m"
    lw  = entry.D == 0.45 ? 2.5 : 1.5
    plot!(p_kmoor_ref, entry.factors, entry.KT_kmoor,
        label=lbl, lw=lw, color=entry.col)
end

vline!(p_kmoor_ref, [1.151], color=:gray, ls=:dot, lw=0.8, label="k/Khs11 = 1.151")
display(p_kmoor_ref)
savefig(p_kmoor_ref, joinpath(plot_dir, "case3_mooring_stiffness_refined.png"))

# ─────────────────────────────────────────────────────────────────────────────
# 9. MOORED DIAMETER SWEEP — using refined optimal stiffness per diameter
# ─────────────────────────────────────────────────────────────────────────────
moored_results_refined = []

@time for (ue, ke) in zip(unmoored_results, kmoor_refined_results)
    res_m = run_frequency_sweep(k_sweep, ue.params, ue.sp, ue.reg, ue.meas;
        params_override=(
            Krb = ComplexF64[ke.k_moor_opt 0.0; 0.0 0.0],
            Crb = ComplexF64[0.0           0.0; 0.0 0.0]
        ),
        label="D=$(ue.D)m moored refined")

    push!(moored_results_refined, (D=ue.D, col=ue.col, res=res_m,
                                   k_opt_factor=ke.k_moor_opt/ke.Khs11))
end

# ── Plot refined moored results ───────────────────────────────────────────────
p_KR_ref = plot(xlabel="ω [rad/s]", ylabel=L"K_{R}", 
                title="Reflection", ylims=(0,1.1), legend=:bottomleft)
p_KT_ref = plot(xlabel="ω [rad/s]", ylabel=L"K_{T}", 
                title="Transmission", ylims=(0,1.1))
p_rao_ref = plot(xlabel="ω [rad/s]", ylabel="RAO heave", 
                 title="Heave RAO")

vline!(p_KR_ref,  [ωp], color=:gray, ls=:dash, lw=0.8, label="ωp")
vline!(p_KT_ref,  [ωp], color=:gray, ls=:dash, lw=0.8, label="ωp")
vline!(p_rao_ref, [ωp], color=:gray, ls=:dash, lw=0.8, label="ωp")

for entry in moored_results_refined
    lbl = "D = $(entry.D) m" * (entry.D == 0.45 ? " ★" : "")
    lw  = entry.D == 0.45 ? 2.5 : 1.5

    plot!(p_KR_ref,  entry.res.ω, entry.res.K_R,       label=lbl, lw=lw, color=entry.col)
    plot!(p_KT_ref,  entry.res.ω, entry.res.K_T,       label=lbl, lw=lw, color=entry.col)
    plot!(p_rao_ref, entry.res.ω, entry.res.RAO_heave, label=lbl, lw=lw, color=entry.col)

    i_min = argmin(entry.res.K_T)
    i_p   = argmin(abs.(entry.res.ω .- ωp))
    @printf("D = %.2f m  k_opt/Khs11 = %.4f  K_T(ωp) = %.4f  K_R(ωp) = %.4f  ω_min = %.3f\n",
        entry.D, entry.k_opt_factor, entry.res.K_T[i_p], entry.res.K_R[i_p], entry.res.ω[i_min])
end

p_D_panel_ref = plot(p_KR_ref, p_KT_ref, p_rao_ref,
    layout=(1,3), size=(1200,400),
    # plot_title="Case 3 — diameter sweep (moored, refined)",
    bottom_margin=7Plots.mm,
    left_margin=7Plots.mm,
    right_margin=7Plots.mm)
display(p_D_panel_ref)
savefig(p_D_panel_ref, joinpath(plot_dir, "case3_D_sweep_moored_refined.png"))

# ── Save all sweep results ────────────────────────────────────────────────────

jldopen(joinpath(jld2_dir_case3, "case3_results.jld2"), "w") do f

    # Unmoored diameter sweep
    for entry in unmoored_results
        d_str = replace(string(entry.D), "." => "p")
        g = JLD2.Group(f, "unmoored/D$(d_str)")
        g["D"]          = entry.D
        g["Khs11"]      = entry.Khs11
        g["m_eff"]      = entry.m_eff
        g["ωn"]         = entry.ωn
        g["k_moor_opt"] = entry.k_moor_opt
        g["ω"]          = entry.res.ω
        g["K_R"]        = entry.res.K_R
        g["K_T"]        = entry.res.K_T
        g["PErr"]       = entry.res.PErr
        g["RAO_heave"]  = entry.res.RAO_heave
    end

    # Coarse mooring stiffness sweep
    for entry in kmoor_sweep_results
        d_str = replace(string(entry.D), "." => "p")
        g = JLD2.Group(f, "kmoor/D$(d_str)")
        g["D"]          = entry.D
        g["Khs11"]      = entry.Khs11
        g["k_moor_opt"] = entry.k_moor_opt
        g["factors"]    = entry.factors
        g["KT_kmoor"]   = entry.KT_kmoor
        g["KR_kmoor"]   = entry.KR_kmoor
    end

    # Refined mooring stiffness sweep
    for entry in kmoor_refined_results
        d_str = replace(string(entry.D), "." => "p")
        g = JLD2.Group(f, "kmoor_refined/D$(d_str)")
        g["D"]          = entry.D
        g["Khs11"]      = entry.Khs11
        g["k_moor_opt"] = entry.k_moor_opt
        g["factors"]    = entry.factors
        g["KT_kmoor"]   = entry.KT_kmoor
        g["KR_kmoor"]   = entry.KR_kmoor
    end

    # Moored diameter sweep (refined optimal stiffness)
    for entry in moored_results_refined
        d_str = replace(string(entry.D), "." => "p")
        g = JLD2.Group(f, "moored_refined/D$(d_str)")
        g["D"]            = entry.D
        g["k_opt_factor"] = entry.k_opt_factor
        g["ω"]            = entry.res.ω
        g["K_R"]          = entry.res.K_R
        g["K_T"]          = entry.res.K_T
        g["PErr"]         = entry.res.PErr
        g["RAO_heave"]    = entry.res.RAO_heave
    end
end

println("Case 3 results saved to: ", joinpath(jld2_dir_case3, "case3_results.jld2"))

# ─────────────────────────────────────────────────────────────────────────────
# 10. WILLEMSPOLDER MOORING STIFFNESS CALCULATION
# Derives k_heave from rope geometry 
# ─────────────────────────────────────────────────────────────────────────────

# Rope data from mooring design (Table 12, ap-w1 to ap-w7)
rope_lengths   = [48.11, 54.70, 50.49, 57.50, 35.78, 46.24, 18.65]  # m
angles_v_min   = [13.3,  14.9,  11.6,  15.8,   8.8,  14.0,   1.4]  # deg
angles_v_max   = [24.0,  24.3,  21.9,  24.7,  23.4,  25.1,  29.7]  # deg

# Axial stiffness from Gleistein X-Permanent Flex 28mm datasheet
# At 10% of break load (307.5 kN): elongation = 0.3%
F_10pct  = 0.1 * 307.5e3   # N
eps_10   = 0.003            # -
EA       = F_10pct / eps_10 # N

# Heave stiffness per anchor: k_vert = (EA/L) * sin²(mean_angle)
k_vert = [EA / L * sind((amin + amax) / 2)^2
          for (L, amin, amax) in zip(rope_lengths, angles_v_min, angles_v_max)]

k_provided_total = sum(k_vert)  # N/m (total for 150m pipe)
L_pipe = 150.0                  # m
k_heave_willemspolder = k_provided_total / L_pipe  # N/m per unit length

# Required optimal stiffness for reference.
# K_hs11_per_m is the per-unit-length hydrostatic stiffness Khs11 for
# D = 0.45 m, read from the unmoored diameter sweep above (case3_results.jld2,
# unmoored/D0p45/Khs11) rather than recomputed here.
K_hs11_per_m  = 4162.0   # N/m/m, D = 0.45 m
k_opt_per_m   = 1.129 * K_hs11_per_m   # 1.129 = optimal k/Khs11 ratio for D = 0.45 m, from the refined mooring sweep (tab:case3_kmoor_refined)

@printf("\n=== Willemspolder mooring stiffness ===\n")
@printf("EA (Gleistein 28mm):         %.0f kN\n",   EA/1e3)
@printf("Total k_provided:            %.1f kN/m\n", k_provided_total/1e3)
@printf("k_heave per unit length:     %.1f N/m/m\n", k_heave_willemspolder)
@printf("k_opt per unit length:       %.1f N/m/m\n", k_opt_per_m)
@printf("Ratio provided/required:     %.2f\n",       k_heave_willemspolder/k_opt_per_m)
@printf("Shortfall factor:            %.1f x\n",     k_opt_per_m/k_heave_willemspolder)

# Rebuild D=0.45m FE objects only
msh_path_45 = "meshes/case3_diametersweep/case3_D_45.msh"
model_45    = Thesis_FPVBarrier.build_model(msh_path_45)
params_45   = Thesis_FPVBarrier.build_case3_params(model_45)
reg_45      = Thesis_FPVBarrier.build_case3_regions(model_45, params_45)
meas_45     = Thesis_FPVBarrier.build_case3_measures(reg_45; order=2)
sp_45       = Thesis_FPVBarrier.build_case3_spaces(reg_45, params_45; order=2)

res_willemspolder = run_frequency_sweep(k_sweep, params_45, sp_45, reg_45, meas_45;
    params_override=(
        Krb = ComplexF64[k_heave_willemspolder 0.0; 0.0 0.0],
        Crb = ComplexF64[0.0                   0.0; 0.0 0.0]
    ),
    label="D=0.45m Willemspolder as-built mooring")

i_p = argmin(abs.(res_willemspolder.ω .- ωp))
@printf("K_T(ωp) = %.4f  K_R(ωp) = %.4f  ΔE_L = %.2f%%  PErr = %.3f%%\n",
    res_willemspolder.K_T[i_p], res_willemspolder.K_R[i_p],
    (1 - res_willemspolder.K_T[i_p]^2)*100, res_willemspolder.PErr[i_p])

jldopen(joinpath(jld2_dir_case3, "case3_results.jld2"), "a+") do f
        # Delete existing group if present
    if haskey(f, "willemspolder_mooring/D0p45")
        delete!(f, "willemspolder_mooring/D0p45")
    end
    g = JLD2.Group(f, "willemspolder_mooring/D0p45")
    g["ω"]               = res_willemspolder.ω
    g["K_T"]             = res_willemspolder.K_T
    g["K_R"]             = res_willemspolder.K_R
    g["PErr"]            = res_willemspolder.PErr
    g["RAO_heave"]       = res_willemspolder.RAO_heave
    g["k_heave_per_m"]   = k_heave_willemspolder
end

# Load unmoored and optimal moored for comparison
data3 = load(joinpath(jld2_dir_case3, "case3_results.jld2"))

ω_u   = data3["unmoored/D0p45/ω"]
KT_u  = abs.(data3["unmoored/D0p45/K_T"])
KR_u  = abs.(data3["unmoored/D0p45/K_R"])
RAO_u = data3["unmoored/D0p45/RAO_heave"]

ω_m   = data3["moored_refined/D0p45/ω"]
KT_m  = abs.(data3["moored_refined/D0p45/K_T"])
KR_m  = abs.(data3["moored_refined/D0p45/K_R"])
RAO_m = data3["moored_refined/D0p45/RAO_heave"]

ω_w   = res_willemspolder.ω
KT_w  = abs.(res_willemspolder.K_T)
KR_w  = abs.(res_willemspolder.K_R)
RAO_w = res_willemspolder.RAO_heave

# Plot
p_KR = plot(xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel=L"K_R", ylims=(0,1.1), title="Reflection", legend=:left)
p_KT = plot(xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel=L"K_T", ylims=(0,1.1), title="Transmission")
p_RAO = plot(xlabel=L"\omega \ \mathrm{[rad/s]}", ylabel="RAO heave", title="Heave RAO")

for (ω, KR, KT, RAO, lbl, col, lw) in [
    (ω_u, KR_u, KT_u, RAO_u, "Unmoored",                :royalblue, 1.5),
    (ω_w, KR_w, KT_w, RAO_w, "Willemspolder as-built",  :orange,    2.0),
    (ω_m, KR_m, KT_m, RAO_m, "Optimal mooring",          :red,       2.0),
]
    plot!(p_KR,  ω, KR,  label=lbl, lw=lw, color=col)
    plot!(p_KT,  ω, KT,  label=lbl, lw=lw, color=col)
    plot!(p_RAO, ω, RAO, label=lbl, lw=lw, color=col)
end

vline!(p_KR,  [ωp], color=:gray, ls=:dash, lw=0.8, label=L"\omega_p")
vline!(p_KT,  [ωp], color=:gray, ls=:dash, lw=0.8, label=L"\omega_p")
vline!(p_RAO, [ωp], color=:gray, ls=:dash, lw=0.8, label=L"\omega_p")

p_panel = plot(p_KR, p_KT, p_RAO,
    layout=(1,3), size=(1200,400),
    bottom_margin=5Plots.mm,
    left_margin=5Plots.mm,
    right_margin=5Plots.mm,
    top_margin=5Plots.mm)

display(p_panel)
savefig(p_panel, joinpath(plot_dir, "case3_willemspolder_mooring_comparison.png"))


# ── Load all sweep results ────────────────────────────────────────────────────
# Use this block instead of rerunning sweeps for plotting or analysis after 
# restarting the Julia session

function load_case3_results(jld2_path)
    f = load(jld2_path)
    D_vals = [0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

    unmoored_loaded = []
    for D in D_vals
        d_str = replace(string(D), "." => "p")
        g = f["unmoored/D$(d_str)"]
        push!(unmoored_loaded, (
            D       = g["D"],
            Khs11   = g["Khs11"],
            m_eff   = g["m_eff"],
            ωn      = g["ωn"],
            k_moor_opt = g["k_moor_opt"],
            res     = (ω=g["ω"], K_R=g["K_R"], K_T=g["K_T"],
                       PErr=g["PErr"], RAO_heave=g["RAO_heave"])
        ))
    end

    kmoor_loaded = []
    for D in D_vals
        d_str = replace(string(D), "." => "p")
        g = f["kmoor/D$(d_str)"]
        push!(kmoor_loaded, (
            D           = g["D"],
            Khs11       = g["Khs11"],
            k_moor_opt  = g["k_moor_opt"],
            factors     = g["factors"],
            KT_kmoor    = g["KT_kmoor"],
            KR_kmoor    = g["KR_kmoor"]
        ))
    end

    kmoor_refined_loaded = []
    for D in D_vals
        d_str = replace(string(D), "." => "p")
        g = f["kmoor_refined/D$(d_str)"]
        push!(kmoor_refined_loaded, (
            D           = g["D"],
            k_factors   = g["k_factors"],
            K_T_vals    = g["K_T_vals"],
            K_R_vals    = g["K_R_vals"],
            k_moor_opt  = g["k_moor_opt"],
            Khs11       = g["Khs11"]
        ))
    end

    moored_loaded = []
    for D in D_vals
        d_str = replace(string(D), "." => "p")
        g = f["moored_refined/D$(d_str)"]
        push!(moored_loaded, (
            D            = g["D"],
            k_opt_factor = g["k_opt_factor"],
            res          = (ω=g["ω"], K_R=g["K_R"], K_T=g["K_T"],
                            PErr=g["PErr"], RAO_heave=g["RAO_heave"])
        ))
    end

    return (unmoored=unmoored_loaded, kmoor=kmoor_loaded,
            kmoor_refined=kmoor_refined_loaded, moored=moored_loaded)
end

# case3_data = load_case3_results(joinpath(jld2_dir_case3, "case3_results.jld2"))
# unmoored_results        = case3_data.unmoored
# kmoor_results           = case3_data.kmoor
# kmoor_refined_results   = case3_data.kmoor_refined
# moored_results_refined  = case3_data.moored

println("\nAll Case 3 postprocessing complete. Plots saved to: ", plot_dir)