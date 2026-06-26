# spectral_weighting.jl
#
# Computes the spectrally-weighted load reduction ΔE_L,spec:
# the JONSWAP-weighted average of the yield indicator 1 - K_T²(ω) over
# the Willemspolder sea state, rather than evaluating it at ωp alone.
# Loads the saved frequency sweeps from case1_postprocess.jl and
# case3_postprocess.jl (results/jld2/case1, results/jld2/case3) and
# compares Case 1 (Lb = 4.6 m) against Case 3 unmoored and optimally
# moored (D = 0.45 m).

using JLD2
using Interpolations
using Plots
using LaTeXStrings
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────────────────────────
base_dir  = joinpath(@__DIR__, "..")
plot_dir  = joinpath(base_dir, "results/plots/energy_yield")
jld2_dir1 = joinpath(base_dir, "results/jld2/case1")
jld2_dir3 = joinpath(base_dir, "results/jld2/case3")
mkpath(plot_dir)

default(fontfamily="DejaVu Sans", guidefontsize=11, tickfontsize=9,
        legendfontsize=8, titlefontsize=14)

# ─────────────────────────────────────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
# JONSWAP spectral shape parameters for the Willemspolder sea state
ωp     = 2.581
γ      = 3.3
α_j    = 0.0081
β_j    = 1.25
ω_low  = 0.5 * ωp
ω_high = 2.0 * ωp
ω_fine = range(ω_low, ω_high, length=1000)

# ─────────────────────────────────────────────────────────────────────────────
# JONSWAP SPECTRUM
# ─────────────────────────────────────────────────────────────────────────────
function jonswap(ω, ωp, γ, α, β_j)
    σ = ω <= ωp ? 0.07 : 0.09
    r = exp(-(ω - ωp)^2 / (2 * σ^2 * ωp^2))
    return α * 9.81^2 / ω^5 * exp(-β_j * (ωp/ω)^4) * γ^r
end

S_fine = [jonswap(ω, ωp, γ, α_j, β_j) for ω in ω_fine]
S_norm = S_fine / maximum(S_fine)

# ─────────────────────────────────────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────
data  = load(joinpath(jld2_dir1, "case1_results.jld2"))
data2 = load(joinpath(jld2_dir3, "case3_results.jld2"))

ω_c1   = data["Lb/4.6/ω"]
KT_c1  = abs.(data["Lb/4.6/K_T"])

ω_c3u  = data2["unmoored/D0p45/ω"]
KT_c3u = abs.(data2["unmoored/D0p45/K_T"])

ω_c3m  = data2["moored_refined/D0p45/ω"]
KT_c3m = abs.(data2["moored_refined/D0p45/K_T"])

# ─────────────────────────────────────────────────────────────────────────────
# INTERPOLATE KT ONTO FINE GRID
# ─────────────────────────────────────────────────────────────────────────────
function interp_KT(ω_vec, KT_vec, ω_fine, ω_low, ω_high)
    idx = ω_low .<= ω_vec .<= ω_high
    itp = LinearInterpolation(ω_vec[idx], KT_vec[idx], extrapolation_bc=Flat())
    return [itp(ω) for ω in ω_fine]
end

IY_c1  = 1 .- interp_KT(ω_c1,  KT_c1,  ω_fine, ω_low, ω_high).^2
IY_c3u = 1 .- interp_KT(ω_c3u, KT_c3u, ω_fine, ω_low, ω_high).^2
IY_c3m = 1 .- interp_KT(ω_c3m, KT_c3m, ω_fine, ω_low, ω_high).^2

# ─────────────────────────────────────────────────────────────────────────────
# SPECTRAL WEIGHTING
# ─────────────────────────────────────────────────────────────────────────────
function delta_EL(IY, S, ω_fine)
    dω = step(ω_fine)
    return sum(S .* IY) * dω / (sum(S) * dω)
end

ΔEL_c1  = delta_EL(IY_c1,  S_fine, ω_fine)
ΔEL_c3u = delta_EL(IY_c3u, S_fine, ω_fine)
ΔEL_c3m = delta_EL(IY_c3m, S_fine, ω_fine)

@printf("ΔE_L,spec Case 1 (Lb=4.6m)          = %.2f%%\n", ΔEL_c1*100)
@printf("ΔE_L,spec Case 3 unmoored (D=0.45m)  = %.2f%%\n", ΔEL_c3u*100)
@printf("ΔE_L,spec Case 3 moored   (D=0.45m)  = %.2f%%\n", ΔEL_c3m*100)

# ─────────────────────────────────────────────────────────────────────────────
# PLOT
# ─────────────────────────────────────────────────────────────────────────────
p = plot(
    xlabel  = L"\omega \ \mathrm{[rad/s]}",
    ylabel  = L"I_\mathrm{yield}(\omega) = 1 - K_T^2(\omega) \ \mathrm{[-]} \quad / \quad S(\omega)/S_\mathrm{max} \ \mathrm{[-]}",
    legend  = :topright,
    ylims   = (0, 1.05),
    xlims   = (ω_low, ω_high),
    yticks  = 0:0.2:1.0,
    grid    = true,
    size    = (800, 420),
    top_margin    = 7Plots.mm,
    bottom_margin = 7Plots.mm,
    left_margin   = 7Plots.mm,
    right_margin  = 7Plots.mm,
)

# JONSWAP shaded background
plot!(p, collect(ω_fine), S_norm,
    fillrange = 0,
    fillalpha = 0.15,
    fillcolor = :steelblue,
    linecolor = :steelblue,
    linestyle = :dash,
    linewidth = 1.2,
    label     = L"S(\omega)/S_\mathrm{max}",
)

# yield indicator curves
plot!(p, collect(ω_fine), IY_c1,
    linewidth = 2,
    color     = :royalblue,
    label     = L"Case 1, $L_b = 4.6$ m " * @sprintf("(ΔEL = %.1f%%)", ΔEL_c1*100),
)

plot!(p, collect(ω_fine), IY_c3u,
    linewidth = 2,
    color     = :red,
    label     = L"Case 3 unmoored, $D = 0.45$ m " * @sprintf("(ΔEL = %.1f%%)", ΔEL_c3u*100),
)

plot!(p, collect(ω_fine), IY_c3m,
    linewidth = 2,
    color     = :green,
    label     = L"Case 3 moored opt, $D = 0.45$ m " * @sprintf("(ΔEL = %.1f%%)", ΔEL_c3m*100),
)

# peak frequency
vline!(p, [ωp],
    color     = :black,
    linestyle = :dash,
    linewidth = 0.8,
    label     = L"\omega_p = 2.581 \ \mathrm{rad/s}",
)

display(p)
savefig(p, joinpath(plot_dir, "yield_indicator_spectral.png"))
println("Figure saved to: ", joinpath(plot_dir, "yield_indicator_spectral.png"))