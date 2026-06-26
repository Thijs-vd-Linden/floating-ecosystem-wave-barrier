# case2_time_postprocess.jl
#
# Post-processing for Case 2: reloads the αp sweep results saved by
# case2_sweep_L4.jl and case2_sweep_L18.jl (one JLD2 file per αp value,
# per beam length) and produces every Case 2 figure used in the
# thesis:
#
#   - Combined free-surface + beam-deflection snapshot at t = 17T,
#     across the full αp sweep (one figure per beam length)
#   - Max beam deflection over time, across the full αp sweep (one
#     figure per beam length)
#   - Beam deflection spatial profile at t = 17T, across the full αp
#     sweep (one figure per beam length)
#
# KT/KR extraction was attempted via Mansard-Funke decomposition but is
# not used in the thesis: outlet-boundary reflections corrupt the
# downstream probe signal with a standing wave, making the
# decomposition unreliable. Results are reported qualitatively
# instead, via the figures above.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "wave_extraction.jl"))
using Plots
using JLD2
using LaTeXStrings

# ── Plot defaults ─────────────────────────────────────────────────────────────
default(fontfamily="DejaVu Sans", guidefontsize=11, tickfontsize=9,
        legendfontsize=8, titlefontsize=14)

# ── Settings ──────────────────────────────────────────────────────────────────
jld2_dir_L4  = "results/jld2/case2/L4"
jld2_dir_L18 = "results/jld2/case2/L18"
plot_dir_L4  = "results/plots/case2_timedomain/postprocess/L4"
plot_dir_L18 = "results/plots/case2_timedomain/postprocess/L18"
mkpath(plot_dir_L4)
mkpath(plot_dir_L18)

αp_sweep  = [0.00, 0.01, 0.1, 0.5, 1.0, 2.5, 5.0]
snap_T    = 17.0   # snapshot time [T]

# Sequential colormap for αp progression
colors = cgrad(:viridis, length(αp_sweep), categorical=true)

# ── Helper: string for αp value ───────────────────────────────────────────────
function ap_to_str(αp)
    replace(string(round(αp, digits=2)), "." => "p")
end

# ── Helper: load JLD2 and return named tuple ──────────────────────────────────
function load_jld2(path)
    data = load(path)
    return (
        κ_profiles = data["κ_profiles"],
        xs_κ       = data["xs_κ"],
        xs_κ_raw   = data["xs_κ_raw"],
        ts         = data["ts"],
        αp         = data["αp"],
        LΩ         = data["LΩ"],
        xb₀        = data["xb₀"],
        xb₁        = data["xb₁"],
        xdᵢₙ       = data["xdᵢₙ"],
        xdₒᵤₜ      = data["xdₒᵤₜ"],
        η₀         = data["η₀"],
        T          = data["T"],
        ω          = data["ω"],
        k          = data["k"],
        κ_T1_t     = data["κ_T1_t"],
        κ_T2_t     = data["κ_T2_t"],
        κ_T3_t     = data["κ_T3_t"],
        κ_R1_t     = data["κ_R1_t"],
        κ_R2_t     = data["κ_R2_t"],
        κ_R3_t     = data["κ_R3_t"],
        κ_max_t    = data["κ_max_t"],
        η_max_t    = data["η_max_t"],
        x_T1       = data["x_T1"],
        x_T2       = data["x_T2"],
        x_T3       = data["x_T3"],
        x_R1       = data["x_R1"],
        x_R2       = data["x_R2"],
        x_R3       = data["x_R3"],
        i_T1       = data["i_T1"],
        i_T2       = data["i_T2"],
        i_T3       = data["i_T3"],
        i_R1       = data["i_R1"],
        i_R2       = data["i_R2"],
        i_R3       = data["i_R3"],
        η_profiles = haskey(data, "η_profiles") ? data["η_profiles"] : nothing,
        xs_η       = haskey(data, "xs_η")       ? data["xs_η"]       : nothing,
    )
end

# ── Combined snapshot figure — Lb = 4.6 m ────────────────────────────────────
println("Generating Lb = 4.6 m snapshot figure...")
plt_L4 = plot(
    xlabel        = L"x / L_{\Omega}",
    ylabel        = L"\cdot / \eta_0",
    title         = "Free surface and beam profiles at \$t = $(Int(snap_T))T\$ — \$L_b = 4.6\$ m",
    framestyle    = :box, grid = true,
    size          = (1000, 400),
    top_margin    = 2Plots.mm,
    bottom_margin =5Plots.mm,
    left_margin   =5Plots.mm,
    right_margin  =5Plots.mm,
    legend        =:outerright
)

d_last_L4 = nothing
for (i, αp) in enumerate(αp_sweep)
    ap_str    = ap_to_str(αp)
    jld2_path = joinpath(jld2_dir_L4, "case2_L4_ap$(ap_str).jld2")
    if !isfile(jld2_path)
        println("  Missing: ", jld2_path, " — skipping")
        continue
    end
    d = load_jld2(jld2_path)
    d_last_L4 = d
    i_snap = argmin(abs.(d.ts .- snap_T*d.T))
    κ_snap = (d.κ_profiles[i_snap])

    # κ left of beam
    mask_L = d.xs_κ .< d.xb₀
    plot!(plt_L4, d.xs_κ[mask_L] ./ d.LΩ, κ_snap[mask_L],
        label = L"\alpha_p = %$(αp) \, \mathrm{s}^{-1}",
        lw    = 1.5,
        color = colors[i]
    )
    # κ right of beam — no label (same color continues)
    mask_R = d.xs_κ .> d.xb₁
    plot!(plt_L4, d.xs_κ[mask_R] ./ d.LΩ, κ_snap[mask_R],
        label = false,
        lw    = 1.5,
        color = colors[i]
    )

    # η beam deflection overlay
    if d.η_profiles !== nothing && d.xs_η !== nothing
        plot!(plt_L4, d.xs_η ./ d.LΩ, (d.η_profiles[i_snap]),
            label = false,
            lw    = 1.5,
            ls    = :dash,
            color = colors[i]
        )
    end
end

if d_last_L4 !== nothing
    d = d_last_L4
    vline!(plt_L4, [d.xb₀ / d.LΩ], ls=:dot, color=:blue, lw=1.0, label=L"x_{b,0}")
    vline!(plt_L4, [d.xb₁ / d.LΩ], ls=:dot, color=:blue, lw=1.0, label=false)
    hline!(plt_L4, [1.0], ls=:dash, color=:black, lw=1.0, label=L"\eta_0")
    plot!(plt_L4, [NaN], [NaN], lw=1.5, ls=:solid, color=:gray, label=L"\kappa")
    plot!(plt_L4, [NaN], [NaN], lw=1.5, ls=:dash,  color=:gray, label=L"\eta")
end

savefig(plt_L4, joinpath(plot_dir_L4, "case2_alpha_sweep_L4.png"))
display(plt_L4)
println("  Saved: case2_alpha_sweep_L4.png")

# ── Combined snapshot figure — Lb = 18.4 m ───────────────────────────────────
println("Generating Lb = 18.4 m snapshot figure...")
plt_L18 = plot(
    xlabel        = L"x / L_{\Omega}",
    ylabel        = L"\cdot / \eta_0",
    title         = "Free surface and beam profiles at \$t = $(Int(snap_T))T\$ — \$L_b = 18.4\$ m",
    framestyle    = :box, grid = true,
    size          = (1000, 400),
    top_margin    = 2Plots.mm,
    bottom_margin =5Plots.mm,
    left_margin   =5Plots.mm,
    right_margin  =5Plots.mm,
    legend        =:outerright
)

d_last_L18 = nothing
for (i, αp) in enumerate(αp_sweep)
    ap_str    = ap_to_str(αp)
    jld2_path = joinpath(jld2_dir_L18, "case2_L18_ap$(ap_str).jld2")
    if !isfile(jld2_path)
        println("  Missing: ", jld2_path, " — skipping")
        continue
    end
    d = load_jld2(jld2_path)
    d_last_L18 = d
    i_snap = argmin(abs.(d.ts .- snap_T*d.T))
    κ_snap = (d.κ_profiles[i_snap])

    # κ left of beam
    mask_L = d.xs_κ .< d.xb₀
    plot!(plt_L18, d.xs_κ[mask_L] ./ d.LΩ, κ_snap[mask_L],
        label = L"\alpha_p = %$(αp) \, \mathrm{s}^{-1}",
        lw    = 1.5,
        color = colors[i]
    )
    # κ right of beam — no label (same color continues)
    mask_R = d.xs_κ .> d.xb₁
    plot!(plt_L18, d.xs_κ[mask_R] ./ d.LΩ, κ_snap[mask_R],
        label = false,
        lw    = 1.5,
        color = colors[i]
    )

    # η beam deflection overlay
    if d.η_profiles !== nothing && d.xs_η !== nothing
        plot!(plt_L18, d.xs_η ./ d.LΩ, (d.η_profiles[i_snap]),
            label = false,
            lw    = 1.5,
            ls    = :dash,
            color = colors[i]
        )
    end
end

if d_last_L18 !== nothing
    d = d_last_L18
    vline!(plt_L18, [d.xb₀ / d.LΩ], ls=:dot, color=:blue, lw=1.0, label=L"x_{b,0}")
    vline!(plt_L18, [d.xb₁ / d.LΩ], ls=:dot, color=:blue, lw=1.0, label=false)
    hline!(plt_L18, [1.0], ls=:dash, color=:black, lw=1.0, label=L"\eta_0")
    plot!(plt_L18, [NaN], [NaN], lw=1.5, ls=:solid, color=:gray, label=L"\kappa")
    plot!(plt_L18, [NaN], [NaN], lw=1.5, ls=:dash,  color=:gray, label=L"\eta")
end

savefig(plt_L18, joinpath(plot_dir_L18, "case2_alpha_sweep_L18.png"))
display(plt_L18)
println("  Saved: case2_alpha_sweep_L18.png")

# ── Beam deflection sweep — Lb = 4.6 m ───────────────────────────────────────
println("Generating Lb = 4.6 m beam deflection figure...")
plt_beam_L4 = plot(
    xlabel     = L"t / T",
    ylabel     = L"\max |\eta| / \eta_0",
    title      = L"Max beam deflection — $L_b = 4.6$ m",
    framestyle = :box, grid = true,
    size       = (800, 400)
)

for (i, αp) in enumerate(αp_sweep)
    ap_str    = ap_to_str(αp)
    jld2_path = joinpath(jld2_dir_L4, "case2_L4_ap$(ap_str).jld2")
    if !isfile(jld2_path)
        println("  Missing: ", jld2_path, " — skipping")
        continue
    end
    d = load_jld2(jld2_path)
    plot!(plt_beam_L4, d.ts ./ d.T, d.η_max_t,
        label = L"\alpha_p = %$(αp) \, \mathrm{s}^{-1}",
        lw    = 1.5,
        color = colors[i])
end

vline!(plt_beam_L4, [snap_T], ls=:dot, color=:black, lw=1.0, label="\$t = $(Int(snap_T))T\$")
savefig(plt_beam_L4, joinpath(plot_dir_L4, "case2_beam_sweep_L4.png"))
display(plt_beam_L4)
println("  Saved: case2_beam_sweep_L4.png")

# ── Beam deflection spatial profile — Lb = 4.6 m ─────────────────────────────
println("Generating Lb = 4.6 m beam deflection spatial profile...")

snap_T_beam = 17

plt_eta_spatial_L4 = plot(
    xlabel     = L"x / L_{\Omega}",
    ylabel     = L"\eta / \eta_0",
    title      = "Beam deflection profiles at \$t = $(Int(snap_T_beam))T\$ — \$L_b = 4.6\$ m",
    framestyle = :box, grid = true,
    size       = (800, 400),
    top_margin    = 2Plots.mm,
    bottom_margin = 5Plots.mm,
    left_margin   = 5Plots.mm,
    right_margin  = 5Plots.mm
)

for (i, αp) in enumerate(αp_sweep)
    ap_str    = ap_to_str(αp)
    jld2_path = joinpath(jld2_dir_L4, "case2_L4_ap$(ap_str).jld2")
    if !isfile(jld2_path)
        println("  Missing: ", jld2_path, " — skipping")
        continue
    end
    d = load_jld2(jld2_path)
    if d.η_profiles === nothing || d.xs_η === nothing
        println("  No η_profiles in $(ap_str) — skipping")
        continue
    end

    i_snap = argmin(abs.(d.ts .- snap_T_beam * d.T))
    println("  αp = $(αp): snap index = $(i_snap), t/T = $(round(d.ts[i_snap]/d.T, digits=3))")

    plot!(plt_eta_spatial_L4, d.xs_η ./ d.LΩ, d.η_profiles[i_snap],
        label = L"\alpha_p = %$(αp) \, \mathrm{s}^{-1}",
        lw    = 1.5,
        color = colors[i]
    )
end

hline!(plt_eta_spatial_L4, [0.0], ls=:dash, color=:black, lw=0.8, label=false)

savefig(plt_eta_spatial_L4, joinpath(plot_dir_L4, "case2_eta_spatial_L4.png"))
display(plt_eta_spatial_L4)
println("  Saved: case2_eta_spatial_L4.png")

# ── Beam deflection spatial profile — Lb = 18.4 m ────────────────────────────
println("Generating Lb = 18.4 m beam deflection spatial profile...")

snap_T_beam = 17.0  # can adjust to any half-period for best shape

plt_eta_spatial_L18 = plot(
    xlabel     = L"x / L_{\Omega}",
    ylabel     = L"\eta / \eta_0",
    title      = "Beam deflection profiles at \$t = $(Int(snap_T_beam))T\$ — \$L_b = 18.4\$ m",
    framestyle = :box, grid = true,
    size       = (800, 400),
    top_margin    = 2Plots.mm,
    bottom_margin = 5Plots.mm,
    left_margin   = 5Plots.mm,
    right_margin  = 5Plots.mm
)

for (i, αp) in enumerate(αp_sweep)
    ap_str    = ap_to_str(αp)
    jld2_path = joinpath(jld2_dir_L18, "case2_L18_ap$(ap_str).jld2")
    if !isfile(jld2_path)
        println("  Missing: ", jld2_path, " — skipping")
        continue
    end
    d = load_jld2(jld2_path)
    if d.η_profiles === nothing || d.xs_η === nothing
        println("  No η_profiles in $(ap_str) — skipping")
        continue
    end

    i_snap = argmin(abs.(d.ts .- snap_T_beam * d.T))
    println("  αp = $(αp): snap index = $(i_snap), t/T = $(round(d.ts[i_snap]/d.T, digits=3))")

    plot!(plt_eta_spatial_L18, d.xs_η ./ d.LΩ, d.η_profiles[i_snap],
        label = L"\alpha_p = %$(αp) \, \mathrm{s}^{-1}",
        lw    = 1.5,
        color = colors[i]
    )
end

hline!(plt_eta_spatial_L18, [0.0], ls=:dash, color=:black, lw=0.8, label=false)

savefig(plt_eta_spatial_L18, joinpath(plot_dir_L18, "case2_eta_spatial_L18.png"))
display(plt_eta_spatial_L18)
println("  Saved: case2_eta_spatial_L18.png")

# ── Beam deflection sweep — Lb = 18.4 m ──────────────────────────────────────
println("Generating Lb = 18.4 m beam deflection figure...")
plt_beam_L18 = plot(
    xlabel     = L"t / T",
    ylabel     = L"\max |\eta| / \eta_0",
    title      = L"Max beam deflection — $L_b = 18.4$ m",
    framestyle = :box, grid = true,
    size       = (800, 400)
)

for (i, αp) in enumerate(αp_sweep)
    ap_str    = ap_to_str(αp)
    jld2_path = joinpath(jld2_dir_L18, "case2_L18_ap$(ap_str).jld2")
    if !isfile(jld2_path)
        println("  Missing: ", jld2_path, " — skipping")
        continue
    end
    d = load_jld2(jld2_path)
    plot!(plt_beam_L18, d.ts ./ d.T, d.η_max_t,
        label = L"\alpha_p = %$(αp) \, \mathrm{s}^{-1}",
        lw    = 1.5,
        color = colors[i])
end

vline!(plt_beam_L18, [snap_T], ls=:dot, color=:black, lw=1.0, label="\$t = $(Int(snap_T))T\$")
savefig(plt_beam_L18, joinpath(plot_dir_L18, "case2_beam_sweep_L18.png"))
display(plt_beam_L18)
println("  Saved: case2_beam_sweep_L18.png")

println("\nAll postprocessing complete.")