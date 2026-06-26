# case1_postprocess.jl
#
# Post-processing for Case 1: floating elastic beam, frequency domain.
# Requires case1.jl to have been run first (provides params1, reg1,
# meas1, sp1).
#
# Re-solves at the design frequency and at a range of frequencies
# around it (independently of whatever case1.jl already solved), to
# produce every Case 1 figure and table in the thesis:
#
#   1. Spatial wave amplitude and beam deflection at ωp (single solve,
#      base case)
#   2. Frequency sweep grid, shared by the three parameter sweeps below
#   3. Bending stiffness (EI) sweep: SDR17 vs SDR11
#   4. Beam density (ρb) sweep: {62.5, 100, 200, 300} kg/m³
#   5. Beam length (Lb) sweep: {4.6, 9.2, 18.4} m
#   6. Bending moment envelope panel, one solve per beam length at ωp
#
# All sweep results are saved to a single JLD2 file
# (results/jld2/case1/case1_results.jld2) so they can be reloaded for
# further analysis or replotting without rerunning the sweeps — see
# load_case1_results() near the end of the file.

using Plots
using Printf
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces
using LaTeXStrings
using JLD2

plot_dir = "results/plots/case1"
jld2_dir_case1 = "results/jld2/case1"
mkpath(plot_dir)
mkpath(jld2_dir_case1)

default(fontfamily="DejaVu Sans", guidefontsize=11, tickfontsize=9, legendfontsize=8, titlefontsize=14)

(; g, H, η₀, ρw, xb₀, xb₁, Lb, xdᵢₙ, xdₒᵤₜ, ωp, kp) = params1

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: extract (x, complex_value) from a CellField on a 1D triangulation
# ─────────────────────────────────────────────────────────────────────────────
function extract_kappa_cells_case1(κh, reg)
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

function extract_field_on_boundary(fh, trian)
    x_vals = Float64[]
    f_vals = ComplexF64[]
    coords   = get_cell_coordinates(trian)
    dof_vals = get_cell_dof_values(fh)
    for i in 1:num_cells(trian)
        cell_coords = coords[i]
        xm   = sum(cell_coords[j][1] for j in 1:length(cell_coords)) / length(cell_coords)
        dofs = dof_vals[i]
        push!(x_vals, xm)
        push!(f_vals, sum(dofs) / length(dofs))
    end
    idx = sortperm(x_vals)
    return x_vals[idx], f_vals[idx]
end

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: Newton-Raphson dispersion relation solver
# ─────────────────────────────────────────────────────────────────────────────
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

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: Mansard-Funke 3-probe wave decomposition
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
# CORE: single frequency solve
# ─────────────────────────────────────────────────────────────────────────────
function solve_single_frequency_case1(k, params_base, sp, reg, meas; params_override=NamedTuple())
    (; g, H, η₀, xb₀, xb₁, xdᵢₙ, xdₒᵤₜ) = params_base

    ω = sqrt(g * k * tanh(k * H))
    λ = 2π / k
    T = 2π / ω

    params_i = merge(params_base, (k=k, ω=ω, λ=λ, T=T))
    params_i = merge(params_i, params_override)

    # @show params_i.k params_i.ω params_i.xb₀ params_i.xb₁ params_i.EI params_i.H
    # @show params_i.Lb params_i.Ld params_i.Lf params_i.LΩ params_i.xdᵢₙ params_i.xdₒᵤₜ

    op_i       = Thesis_FPVBarrier.build_case1_operator(sp, reg, meas, 2, params_i)
    _, κh_i, ηh_i = solve(op_i)

    x_cells, κ_cells = extract_kappa_cells_case1(κh_i, reg)

    # Upstream 3-probe
    x_up_end     = xb₀ - (xb₀ - xdᵢₙ) * 0.05
    x_p1         = xdᵢₙ + (xb₀ - xdᵢₙ) * 0.30
    Δx_up        = min(λ / 4, (x_up_end - x_p1) / 2.5)
    A_inc, B_ref = decompose_wave(x_cells, κ_cells, k,
                       x_p1, x_p1 + Δx_up, x_p1 + 2Δx_up)
    K_R = abs(B_ref) / abs(A_inc)

    # Downstream 3-probe
    x_dn_start   = xb₁  + (xdₒᵤₜ - xb₁) * 0.10
    x_dn_end     = xdₒᵤₜ - (xdₒᵤₜ - xb₁) * 0.05
    Δx_dn        = min(λ / 4, (x_dn_end - x_dn_start) / 2.5)
    C_trans, _   = decompose_wave(x_cells, κ_cells, k,
                       x_dn_start, x_dn_start + Δx_dn, x_dn_start + 2Δx_dn)
    K_T = abs(C_trans) / abs(A_inc)

    F_drift_nd = 1.0 + K_R - K_T
    PErr       = (1.0 - K_R^2 - K_T^2) * 100

    return (ω=ω, K_R=K_R, K_T=K_T, F_drift=F_drift_nd, PErr=PErr,
            A_inc=A_inc, B_ref=B_ref, C_trans=C_trans,
            x_p1=x_p1, x_p2=x_p1+Δx_up, x_p3=x_p1+2Δx_up,
            x_dn1=x_dn_start, x_dn2=x_dn_start+Δx_dn, x_dn3=x_dn_start+2Δx_dn,
            Δx_up=Δx_up, λ=λ, κh=κh_i, ηh=ηh_i)
end

# ─────────────────────────────────────────────────────────────────────────────
# CORE: full frequency sweep
# ─────────────────────────────────────────────────────────────────────────────
function run_frequency_sweep_case1(k_range, params_base, sp, reg, meas;
                                    params_override=NamedTuple(), label="")
    N = length(k_range)
    ω_vals    = Float64[]
    K_R_vals  = Float64[]
    K_T_vals  = Float64[]
    F_drift   = Float64[]
    PErr_vals = Float64[]

    println("Frequency sweep: $label ($N frequencies)...")
    flush(stdout)

    for (i, k) in enumerate(k_range)
        r = solve_single_frequency_case1(k, params_base, sp, reg, meas;
                                          params_override=params_override)
        push!(ω_vals,    r.ω)
        push!(K_R_vals,  r.K_R)
        push!(K_T_vals,  r.K_T)
        push!(F_drift,   r.F_drift)
        push!(PErr_vals, r.PErr)

        if mod(i, 10) == 0 || i == 1
            @printf("  [%2d/%d]  k=%.3f  ω=%.3f  K_R=%.3f  K_T=%.3f  PErr=%.1f%%\n",
                i, N, k, r.ω, r.K_R, r.K_T, r.PErr)
            flush(stdout)
        end
    end

    println("  Done: $label")
    return (ω=ω_vals, K_R=K_R_vals, K_T=K_T_vals, F_drift=F_drift, PErr=PErr_vals)
end

# ─────────────────────────────────────────────────────────────────────────────
# 1. SPATIAL PLOT at ωp (single frequency, base params)
# ─────────────────────────────────────────────────────────────────────────────
function plot_kappa_eta_case1(κh, ηh, reg, params)
    (; η₀, LΩ, k, xdᵢₙ, xdₒᵤₜ, xb₀, xb₁) = params

    x_fs, κ_vals  = extract_field_on_boundary(κh, reg.Γfs)
    κ_in_vals      = η₀ .* exp.(im .* k .* x_fs)
    κ_r_vals       = κ_vals .- κ_in_vals
    x_str, η_vals  = extract_field_on_boundary(ηh, reg.Γstr)

    κ_abs    = abs.(κ_vals)    ./ η₀
    κ_r_abs  = abs.(κ_r_vals)  ./ η₀
    κ_in_abs = abs.(κ_in_vals) ./ η₀
    η_abs    = abs.(η_vals)    ./ η₀

    mask_left  = x_fs .< xb₀
    mask_right = x_fs .> xb₁

    p1 = plot(
        x_fs[mask_left] ./ LΩ, κ_abs[mask_left],
        label="Free surface", lw=2, color=:black,
        xlabel=L"x / L_{\Omega}", ylabel=L"Normalised amplitude $|\cdot| / \eta_0$",
        legend=:bottomright
    )
    plot!(p1, x_fs[mask_right] ./ LΩ, κ_abs[mask_right],  label="", lw=2, color=:black)
    plot!(p1, x_fs[mask_left]  ./ LΩ, κ_r_abs[mask_left], label="Reflected wave", lw=1.5, color=:royalblue)
    plot!(p1, x_fs[mask_left]  ./ LΩ, κ_in_abs[mask_left], label="Incident wave", lw=1.0, color=:gray, ls=:dot)
    plot!(p1, x_str ./ LΩ, η_abs, label="Beam", lw=2, color=:red)
    vline!(p1, [xdᵢₙ  / LΩ], color=:green, ls=:dash, lw=0.8, label="Damping edges")
    vline!(p1, [xdₒᵤₜ / LΩ], color=:green, ls=:dash, lw=0.8, label="")
    vspan!(p1, [xb₀/LΩ, xb₁/LΩ], alpha=0.08, color=:red, label="Beam region")

    display(p1)
    savefig(p1, joinpath(plot_dir, "case1_kappa_eta_spatial.png"))

    small_mask_left  = (x_fs .>= xdᵢₙ) .& (x_fs .< xb₀)
    small_mask_right = (x_fs .> xb₁)   .& (x_fs .<= xdₒᵤₜ)

    p2 = plot(
        x_fs[small_mask_left] ./ LΩ, κ_abs[small_mask_left],
        label=L"|\kappa| / \eta_0", lw=2, color=:black,
        xlabel=L"x / L_{\Omega}", ylabel=L"$|\kappa| / \eta_0$",
        legend=:bottomright,
        # title="Surface elevation — Case 1"
    )
    plot!(p2, x_fs[small_mask_right] ./ LΩ, κ_abs[small_mask_right], label="", lw=2, color=:black)
    plot!(p2, x_fs[small_mask_left]  ./ LΩ, κ_r_abs[small_mask_left],  label=L"|\kappa_r| / \eta_0", lw=1.5, color=:royalblue)
    plot!(p2, x_fs[small_mask_left]  ./ LΩ, κ_in_abs[small_mask_left], label=L"|\kappa_{in}| / \eta_0", lw=1.0, color=:gray, ls=:dot)
    plot!(p2, x_str ./ LΩ, η_abs, label=L"|\eta| / \eta_0", lw=2, color=:red)
    vspan!(p2, [xb₀/LΩ, xb₁/LΩ], alpha=0.08, color=:red, label="Beam region")

    display(p2)
    savefig(p2, joinpath(plot_dir, "case1_kappa_eta_spatial_inner.png"))
    return p1
end

r_spatial = solve_single_frequency_case1(kp, params1, sp1, reg1, meas1)
plot_kappa_eta_case1(r_spatial.κh, r_spatial.ηh, reg1, merge(params1, (k=kp, ω=r_spatial.ω, λ=r_spatial.λ, T=2π/r_spatial.ω)))

@printf("  x_p1=%.3f  x_p2=%.3f  x_p3=%.3f  xb₀=%.3f\n", 
    r_spatial.x_p1, r_spatial.x_p2, r_spatial.x_p3, params1.xb₀)

# ─────────────────────────────────────────────────────────────────────────────
# 2. FREQUENCY SWEEP — ω range matching Case 3
# ─────────────────────────────────────────────────────────────────────────────
ω_low  = range(0.5*ωp, 2.0, length=30)
ω_mid  = range(2.0, ωp, length=25)
ω_high = range(ωp, 2.0*ωp, length=40)
ω_sweep = unique(vcat(collect(ω_low), collect(ω_mid), collect(ω_high)))
k_sweep = find_k.(ω_sweep, H)

# ─────────────────────────────────────────────────────────────────────────────
# 3. EI SWEEP — SDR11 and SDR17
# ─────────────────────────────────────────────────────────────────────────────
model_EI  = Thesis_FPVBarrier.build_model("meshes/case1_length_sweep/case1_L_18.msh")
params_EI = Thesis_FPVBarrier.build_case1_params(model_EI)
reg_EI    = Thesis_FPVBarrier.build_case1_regions(model_EI, params_EI)
meas_EI   = Thesis_FPVBarrier.build_case1_measures(reg_EI; order=order)
sp_EI     = Thesis_FPVBarrier.build_case1_spaces(reg_EI; order=order)

@printf("  EI sweep mesh check: Lb=%.2f  xb₀=%.3f  xb₁=%.3f\n", params_EI.Lb, params_EI.xb₀, params_EI.xb₁)

EI_vals  = [params_EI.EI_SDR17, params_EI.EI_SDR11]
EI_names = ["SDR17", "SDR11"]
EI_cols  = [:royalblue, :coral]

EI_results = []

println("\nEI sweep (SDR11 / SDR17) on Lb = 18.4 m mesh...")
@time for (EI_i, name_i, col_i) in zip(EI_vals, EI_names, EI_cols)
    res_i = run_frequency_sweep_case1(k_sweep, params_EI, sp_EI, reg_EI, meas_EI;
        params_override=(EI=EI_i,),
        label="Case 1 $name_i")
    push!(EI_results, (name=name_i, col=col_i, EI=EI_i, res=res_i))
end

# ─────────────────────────────────────────────────────────────────────────────
# 3b. EI SWEEP PLOTS
# ─────────────────────────────────────────────────────────────────────────────
p_KR_EI = plot(xlabel=L"ω [rad/s]", ylabel=L"K_{R}", title="Reflection", ylims=(0,1.1))
p_KT_EI = plot(xlabel=L"ω [rad/s]", ylabel=L"K_{T}", title="Transmission", ylims=(0,1.1))
p_err_EI = plot(xlabel=L"ω [rad/s]", ylabel=L"PErr \ [\%]", title="Energy error")

vline!(p_KR_EI,  [ωp], color=:gray, ls=:dash, lw=0.8, label=L"ω_{p}")
vline!(p_KT_EI,  [ωp], color=:gray, ls=:dash, lw=0.8, label=L"ω_{p}")
vline!(p_err_EI, [ωp], color=:gray, ls=:dash, lw=0.8, label=L"ω_{p}")
hline!(p_err_EI, [0],  color=:black, ls=:dot, lw=0.8, label="")

for entry in EI_results
    plot!(p_KR_EI,  entry.res.ω, entry.res.K_R,  label=entry.name, lw=2, color=entry.col)
    plot!(p_KT_EI,  entry.res.ω, entry.res.K_T,  label=entry.name, lw=2, color=entry.col)
    plot!(p_err_EI, entry.res.ω, entry.res.PErr, label=entry.name, lw=2, color=entry.col)
end

p_EI_panel = plot(p_KR_EI, p_KT_EI, p_err_EI,
    layout=(1,3), size=(1200,350),
    top_margin=7Plots.mm,
    bottom_margin=7Plots.mm,
    left_margin=7Plots.mm,
    right_margin=7Plots.mm)
display(p_EI_panel)
savefig(p_EI_panel, joinpath(plot_dir, "case1_EI_sweep.png"))

# ─────────────────────────────────────────────────────────────────────────────
# 4b. ρb SWEEP — {62.5, 100, 200, 300} kg/m³
# ─────────────────────────────────────────────────────────────────────────────
model_ρb  = Thesis_FPVBarrier.build_model("meshes/case1_length_sweep/case1_L_18.msh")
params_ρb = Thesis_FPVBarrier.build_case1_params(model_ρb)
reg_ρb    = Thesis_FPVBarrier.build_case1_regions(model_ρb, params_ρb)
meas_ρb   = Thesis_FPVBarrier.build_case1_measures(reg_ρb; order=order)
sp_ρb     = Thesis_FPVBarrier.build_case1_spaces(reg_ρb; order=order)

@printf("  ρb sweep mesh check: Lb=%.2f  xb₀=%.3f  xb₁=%.3f\n", params_ρb.Lb, params_ρb.xb₀, params_ρb.xb₁)

ρb_vals  = [62.5, 100.0, 200.0, 300.0]
ρb_names = ["62.5", "100", "200", "300"]
ρb_cols  = [:royalblue, :green, :orange, :red]

ρb_results = []

println("\nρb sweep...")
@time for (ρb_i, name_i, col_i) in zip(ρb_vals, ρb_names, ρb_cols)
    res_i = run_frequency_sweep_case1(k_sweep, params_ρb, sp_ρb, reg_ρb, meas_ρb;
        params_override=(ρb=ρb_i,),
        label="Case 1 ρb=$name_i")
    push!(ρb_results, (name=name_i, col=col_i, ρb=ρb_i, res=res_i))
end

p_KR_ρb  = plot(xlabel=L"ω [rad/s]", ylabel=L"K_{R}", title="Reflection", ylims=(0,1.1), legend=:topleft)
p_KT_ρb  = plot(xlabel=L"ω [rad/s]", ylabel=L"K_{T}", title="Transmission", ylims=(0,1.1), legend=:left)
p_err_ρb = plot(xlabel=L"ω [rad/s]", ylabel=L"PErr \ [\%]", title="Energy error")

vline!(p_KR_ρb,  [ωp], color=:gray, ls=:dash, lw=0.8, label=L"ω_{p}")
vline!(p_KT_ρb,  [ωp], color=:gray, ls=:dash, lw=0.8, label=L"ω_{p}")
vline!(p_err_ρb, [ωp], color=:gray, ls=:dash, lw=0.8, label=L"ω_{p}")
hline!(p_err_ρb, [0],  color=:black, ls=:dot, lw=0.8, label="")

for entry in ρb_results
    lbl = L"\rho_b = %$(entry.name) \ \mathrm{kg/m^3}" * (entry.ρb == 62.5 ? " ★" : "")
    lw  = entry.ρb == 62.5 ? 2.5 : 1.5
    plot!(p_KR_ρb,  entry.res.ω, entry.res.K_R,  label=lbl, lw=lw, color=entry.col)
    plot!(p_KT_ρb,  entry.res.ω, entry.res.K_T,  label=lbl, lw=lw, color=entry.col)
    plot!(p_err_ρb, entry.res.ω, entry.res.PErr, label=lbl, lw=lw, color=entry.col)
end

p_ρb_panel = plot(p_KR_ρb, p_KT_ρb, p_err_ρb,
    layout=(1,3), size=(1200,350),
    top_margin=7Plots.mm, bottom_margin=7Plots.mm, left_margin=7Plots.mm, right_margin=7Plots.mm)
display(p_ρb_panel)
savefig(p_ρb_panel, joinpath(plot_dir, "case1_rho_sweep.png"))

# ─────────────────────────────────────────────────────────────────────────────
# 5. Lb SWEEP — {4.6, 9.2, 18.4} m 
# ─────────────────────────────────────────────────────────────────────────────
Lb_vals   = [4.6, 9.2, 18.4]
Lb_names  = ["4.6", "9.2", "18.4"]
Lb_cols   = [:royalblue, :green, :red]
Lb_meshes = [
    "meshes/case1_length_sweep/case1_L_4.msh",
    "meshes/case1_length_sweep/case1_L_9.msh",
    "meshes/case1_length_sweep/case1_L_18.msh",
]


Lb_results = []

println("\nLb sweep...")
@time for (Lb_i, name_i, col_i, msh_path) in zip(Lb_vals, Lb_names, Lb_cols, Lb_meshes)
    model_i  = Thesis_FPVBarrier.build_model(msh_path)
    params_i = Thesis_FPVBarrier.build_case1_params(model_i)
    @printf("  Lb mesh check: Lb=%.2f  xb₀=%.3f  xb₁=%.3f\n", params_i.Lb, params_i.xb₀, params_i.xb₁)
    reg_i    = Thesis_FPVBarrier.build_case1_regions(model_i, params_i)
    meas_i   = Thesis_FPVBarrier.build_case1_measures(reg_i; order=order)
    sp_i     = Thesis_FPVBarrier.build_case1_spaces(reg_i; order=order)
    res_i    = run_frequency_sweep_case1(k_sweep, params_i, sp_i, reg_i, meas_i;
                   label="Case 1 Lb=$name_i m")
    push!(Lb_results, (name=name_i, col=col_i, Lb=Lb_i, res=res_i))
end

p_KR_Lb  = plot(xlabel=L"ω [rad/s]", ylabel=L"K_{R}", title="Reflection", ylims=(0,1.1), legend=:topleft)
p_KT_Lb  = plot(xlabel=L"ω [rad/s]", ylabel=L"K_{T}", title="Transmission", ylims=(0,1.1), legend=:left)
p_err_Lb = plot(xlabel=L"ω [rad/s]", ylabel=L"PErr \ [\%]", title="Energy error")

vline!(p_KR_Lb,  [ωp], color=:gray, ls=:dash, lw=0.8, label=L"ω_{p}")
vline!(p_KT_Lb,  [ωp], color=:gray, ls=:dash, lw=0.8, label=L"ω_{p}")
vline!(p_err_Lb, [ωp], color=:gray, ls=:dash, lw=0.8, label=L"ω_{p}")
hline!(p_err_Lb, [0],  color=:black, ls=:dot, lw=0.8, label="")

for entry in Lb_results
    lbl = L"L_{b} = %$(entry.name) \ \mathrm{m}" * (entry.Lb == 4.6 ? " ★" : "")
    lw  = entry.Lb == 4.6 ? 2.5 : 1.5
    plot!(p_KR_Lb,  entry.res.ω, entry.res.K_R,  label=lbl, lw=lw, color=entry.col)
    plot!(p_KT_Lb,  entry.res.ω, entry.res.K_T,  label=lbl, lw=lw, color=entry.col)
    plot!(p_err_Lb, entry.res.ω, entry.res.PErr, label=lbl, lw=lw, color=entry.col)
end

p_Lb_panel = plot(p_KR_Lb, p_KT_Lb, p_err_Lb,
    layout=(1,3), size=(1200,350),
    top_margin=7Plots.mm, bottom_margin=7Plots.mm, left_margin=7Plots.mm, right_margin=7Plots.mm)
display(p_Lb_panel)
savefig(p_Lb_panel, joinpath(plot_dir, "case1_Lb_sweep.png"))

println("All plots saved to ", plot_dir)


# ─────────────────────────────────────────────────────────────────────────────
# BENDING MOMENT PANEL — L4, L9, L18 at ωp
# ─────────────────────────────────────────────────────────────────────────────
η₀_s  = 0.249   # significant wave amplitude [m]
EI    = params1.EI  # base params, same EI for all lengths

p_M_list = []

for (msh, lb_name, Lb_i) in zip(Lb_meshes, Lb_names, Lb_vals)
    # Build model at ωp
    model_i  = Thesis_FPVBarrier.build_model(msh)
    params_i = Thesis_FPVBarrier.build_case1_params(model_i)
    reg_i    = Thesis_FPVBarrier.build_case1_regions(model_i, params_i)
    meas_i   = Thesis_FPVBarrier.build_case1_measures(reg_i; order=order)
    sp_i     = Thesis_FPVBarrier.build_case1_spaces(reg_i; order=order)

    # Solve at ωp using params_i (already set to ωp)
    op_i          = Thesis_FPVBarrier.build_case1_operator(sp_i, reg_i, meas_i, order, params_i)
    ϕh_i, κh_i, ηh_i = solve(op_i)

    # Extract η along beam
    x_str, η_vals = extract_field_on_boundary(ηh_i, reg_i.Γstr)

    # Second derivative via central differences
    n     = length(x_str)
    M_vals = Float64[]
    x_mid  = Float64[]
    for j in 2:n-1
        dx   = (x_str[j+1] - x_str[j-1]) / 2
        η_xx = (η_vals[j+1] - 2*η_vals[j] + η_vals[j-1]) / dx^2
        push!(M_vals, abs(EI * η_xx))
        push!(x_mid,  x_str[j])
    end

    # Scale to physical units
    scale  = η₀_s / params_i.η₀
    M_phys = scale .* M_vals ./ 1000  # [kNm/m]

    # Normalise x to [0, 1]
    x_norm = (x_mid .- params_i.xb₀) ./ Lb_i

    # Plot
    p_i = plot(
        x_norm, M_phys,
        xlabel = L"x / L_b",
        ylabel = lb_name == "4.6" ? L"|M| \ [\mathrm{kNm/m}]" : "",
        title  = L"L_b = %$(lb_name) \ \mathrm{m}",
        lw     = 2,
        color  = :black,
        legend = false,
        ylims  = (0, nothing),
        grid   = true,
        framestyle = :box,
    )
    push!(p_M_list, p_i)
end

p_M_panel = plot(p_M_list...,
    layout      = (1, 3),
    size        = (1100, 320),
    top_margin    = 5Plots.mm,
    bottom_margin = 8Plots.mm,
    left_margin   = 8Plots.mm,
    right_margin  = 3Plots.mm,
)

display(p_M_panel)
savefig(p_M_panel, joinpath(plot_dir, "case1_bending_moment_panel.png"))

# ── Save all sweep results ────────────────────────────────────────────────────

jldopen(joinpath(jld2_dir_case1, "case1_results.jld2"), "w") do f

    # k/ω sweep axis
    f["k_sweep"] = k_sweep
    f["ω_sweep"] = ω_sweep

    # EI sweep
    for entry in EI_results
        g = JLD2.Group(f, "EI/$(entry.name)")
        g["EI"]    = entry.EI
        g["ω"]     = entry.res.ω
        g["K_R"]   = entry.res.K_R
        g["K_T"]   = entry.res.K_T
        g["PErr"]  = entry.res.PErr
        g["F_drift"] = entry.res.F_drift
    end

    # ρb sweep
    for entry in ρb_results
        g = JLD2.Group(f, "rhob/$(entry.name)")
        g["rhob"]  = entry.ρb
        g["ω"]     = entry.res.ω
        g["K_R"]   = entry.res.K_R
        g["K_T"]   = entry.res.K_T
        g["PErr"]  = entry.res.PErr
        g["F_drift"] = entry.res.F_drift
    end

    # Lb sweep
    for entry in Lb_results
        g = JLD2.Group(f, "Lb/$(entry.name)")
        g["Lb"]    = entry.Lb
        g["ω"]     = entry.res.ω
        g["K_R"]   = entry.res.K_R
        g["K_T"]   = entry.res.K_T
        g["PErr"]  = entry.res.PErr
        g["F_drift"] = entry.res.F_drift
    end

    # Spatial solve at ωp
    x_fs, κ_vals = extract_field_on_boundary(r_spatial.κh, reg1.Γfs)
    x_str, η_vals = extract_field_on_boundary(r_spatial.ηh, reg1.Γstr)

    f["spatial/x_fs"]   = x_fs
    f["spatial/κ_vals"] = κ_vals        # complex
    f["spatial/x_str"]  = x_str
    f["spatial/η_vals"] = η_vals        # complex
    f["spatial/ω"]      = r_spatial.ω
    f["spatial/K_R"]    = r_spatial.K_R
    f["spatial/K_T"]    = r_spatial.K_T
    f["spatial/x_p1"]   = r_spatial.x_p1
    f["spatial/x_p2"]   = r_spatial.x_p2
    f["spatial/x_p3"]   = r_spatial.x_p3
    f["spatial/x_dn1"]  = r_spatial.x_dn1
    f["spatial/x_dn2"]  = r_spatial.x_dn2
    f["spatial/x_dn3"]  = r_spatial.x_dn3
    f["spatial/λ"]      = r_spatial.λ
end

println("Case 1 results saved to: ", joinpath(jld2_dir_case1, "case1_results.jld2"))


# ── Load all sweep results ────────────────────────────────────────────────────
# Use this block instead of rerunning sweeps for plotting or analysis after 
# restarting the Julia session

function load_case1_results(jld2_path)
    f = load(jld2_path)

    EI_loaded = []
    for name in ["SDR17", "SDR11"]
        g = f["EI/$name"]
        push!(EI_loaded, (
            name = name,
            EI   = g["EI"],
            res  = (ω=g["ω"], K_R=g["K_R"], K_T=g["K_T"],
                    PErr=g["PErr"], F_drift=g["F_drift"])
        ))
    end

    ρb_loaded = []
    for name in ["62.5", "100", "200", "300"]
        g = f["rhob/$name"]
        push!(ρb_loaded, (
            name = name,
            ρb   = g["rhob"],
            res  = (ω=g["ω"], K_R=g["K_R"], K_T=g["K_T"],
                    PErr=g["PErr"], F_drift=g["F_drift"])
        ))
    end

    Lb_loaded = []
    for name in ["4.6", "9.2", "18.4"]
        g = f["Lb/$name"]
        push!(Lb_loaded, (
            name = name,
            Lb   = g["Lb"],
            res  = (ω=g["ω"], K_R=g["K_R"], K_T=g["K_T"],
                    PErr=g["PErr"], F_drift=g["F_drift"])
        ))
    end

    return (EI=EI_loaded, rhob=ρb_loaded, Lb=Lb_loaded,
            k_sweep=f["k_sweep"], ω_sweep=f["ω_sweep"])
end

# ── Uncomment the following lines to load results from the JLD2 file instead of rerunning sweeps ───────────────
# case1_data = load_case1_results(joinpath(jld2_dir_case1, "case1_results.jld2"))
# EI_results  = case1_data.EI
# ρb_results  = case1_data.rhob
# Lb_results  = case1_data.Lb

## ── Print for report ────────────────────────────────────────────────────────────
# Print values at ωp for all sweeps
i_ωp = argmin(abs.(EI_results[1].res.ω .- ωp));
for entry in EI_results
    @printf("EI %s:  ω=%.3f  K_R=%.4f  K_T=%.4f  PErr=%.3f%%\n",
        entry.name, entry.res.ω[i_ωp], entry.res.K_R[i_ωp], entry.res.K_T[i_ωp], entry.res.PErr[i_ωp])
end

i_ωp_ρb = argmin(abs.(ρb_results[1].res.ω .- ωp));
for entry in ρb_results
    @printf("ρb %s:  ω=%.3f  K_R=%.4f  K_T=%.4f  PErr=%.3f%%\n",
        entry.name, entry.res.ω[i_ωp_ρb], entry.res.K_R[i_ωp_ρb], entry.res.K_T[i_ωp_ρb], entry.res.PErr[i_ωp_ρb])
end

i_ωp_Lb = argmin(abs.(Lb_results[1].res.ω .- ωp));
for entry in Lb_results
    @printf("Lb %s:  ω=%.3f  K_R=%.4f  K_T=%.4f  PErr=%.3f%%\n",
        entry.name, entry.res.ω[i_ωp_Lb], entry.res.K_R[i_ωp_Lb], entry.res.K_T[i_ωp_Lb], entry.res.PErr[i_ωp_Lb])
end