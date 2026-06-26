# structural_impl.jl
#
# Computes the spectrally-weighted load reduction ΔE_L,spec, the
# associated load and fatigue damage-rate reduction factors R_load and
# R_fatigue = R_load^n (for the S-N exponents in n_values, per
# EN 1993-1-9), and the resulting reduction in fatigue cycle count
# above the Rayleigh tail, over the design life T_design.
#
# Loads the saved frequency sweeps from case1_postprocess.jl and
# case3_postprocess.jl (results/jld2/case1, results/jld2/case3) and
# compares four configurations: Case 1 (Lb = 4.6 m), Case 3 unmoored,
# Case 3 optimally moored, and Case 3 with the Willemspolder as-built
# mooring stiffness.

using JLD2
using Statistics
using Printf
using LinearAlgebra

# ─────────────────────────────────────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
ω_p      = 2.581       # rad/s, peak frequency
H_s      = 0.498       # m, significant wave height
η0_p     = H_s / 2     # m, peak wave amplitude = Hs/2
γ        = 3.3         # JONSWAP peak enhancement factor
n_values = [3, 5]      # S-N exponents
T_design = 25          # years, design life
α_j      = 0.0081      # JONSWAP spectral energy scale parameter
β_j      = 1.25        # JONSWAP peak shape parameter

# ─────────────────────────────────────────────────────────────────────────────
# JONSWAP SPECTRUM
# ─────────────────────────────────────────────────────────────────────────────
function jonswap(ω, ω_p, γ, α_j, β_j)
    σ = ω <= ω_p ? 0.07 : 0.09
    r = exp(-(ω - ω_p)^2 / (2 * σ^2 * ω_p^2))
    return α_j * 9.81^2 / ω^5 * exp(-β_j * (ω_p/ω)^4) * γ^r
end

# ─────────────────────────────────────────────────────────────────────────────
# LOAD KT DATA
# ─────────────────────────────────────────────────────────────────────────────
jld2_dir1 = "results/jld2/case1"
jld2_dir3 = "results/jld2/case3"

data  = load(joinpath(jld2_dir1, "case1_results.jld2"))
data2 = load(joinpath(jld2_dir3, "case3_results.jld2"))

ω_grid = data["Lb/4.6/ω"]
KT_c1  = abs.(data["Lb/4.6/K_T"])

ω_c3   = data2["unmoored/D0p45/ω"]
KT_un  = clamp.(abs.(data2["unmoored/D0p45/K_T"]),       0.0, 1.0)
KT_mo  = clamp.(abs.(data2["moored_refined/D0p45/K_T"]), 0.0, 1.0)

# willemspolder mooring stiffness comparison data
ω_c3wp  = data2["willemspolder_mooring/D0p45/ω"]
KT_c3wp = clamp.(abs.(data2["willemspolder_mooring/D0p45/K_T"]), 0.0, 1.0)

# Interpolate Case 3 onto Case 1 grid
function interp1(ω_new, ω_old, vals)
    [begin
        i = clamp(searchsortedlast(ω_old, ω), 1, length(ω_old)-1)
        t = (ω - ω_old[i]) / (ω_old[i+1] - ω_old[i])
        vals[i] * (1-t) + vals[i+1] * t
    end for ω in ω_new]
end

KT_un_g = interp1(ω_grid, ω_c3, KT_un)      # Interpolated unmoored Case 3 data onto Case 1 grid
KT_mo_g = interp1(ω_grid, ω_c3, KT_mo)      # Interpolated moored Case 3 data onto Case 1 grid
KT_wp_g = interp1(ω_grid, ω_c3wp, KT_c3wp)  # Interpolated Willemspolder mooring data onto Case 1 grid

# ─────────────────────────────────────────────────────────────────────────────
# JONSWAP ON SWEEP GRID — TRAPEZOIDAL INTEGRATION
# ─────────────────────────────────────────────────────────────────────────────
S_inc  = jonswap.(ω_grid, ω_p, γ, α_j, β_j)
m0_inc = sum(0.5 .* (S_inc[1:end-1] .+ S_inc[2:end]) .* diff(ω_grid))

configs = [
    ("Case 1, Lb=4.6m SDR17",         KT_c1),
    ("Case 3, unmoored D=0.45m",       KT_un_g),
    ("Case 3, moored optimal D=0.45m", KT_mo_g),
    ("Case 3, Willemspolder as-built",      KT_wp_g),
]

# ─────────────────────────────────────────────────────────────────────────────
# SPECTRAL LOAD REDUCTION
# ─────────────────────────────────────────────────────────────────────────────
println("=" ^ 70)
println("SPECTRAL LOAD REDUCTION RESULTS")
println("=" ^ 70)
println(rpad("Configuration", 38),
        rpad("ΔEL_spec [%]", 14),
        rpad("R_load", 10),
        rpad("R_drift", 10))
println("-" ^ 70)

results = []
for (name, KT) in configs
    S_trans  = KT.^2 .* S_inc
    m0_trans = sum(0.5 .* (S_trans[1:end-1] .+ S_trans[2:end]) .* diff(ω_grid))
    ΔEL      = 1 - m0_trans / m0_inc
    R_load   = sqrt(m0_trans / m0_inc)
    R_drift  = m0_trans / m0_inc
    push!(results, (name, ΔEL, R_load, R_drift, KT, S_trans, m0_trans))
    println(rpad(name, 38),
            rpad(round(ΔEL*100,  digits=1), 14),
            rpad(round(R_load,   digits=3), 10),
            rpad(round(R_drift,  digits=3), 10))
end

# ─────────────────────────────────────────────────────────────────────────────
# FATIGUE DAMAGE RATE RATIOS
# ─────────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 70)
println("FATIGUE DAMAGE RATE RATIOS  (D_trans / D_inc = R_load^n)")
println("=" ^ 70)
println(rpad("Configuration", 38), join([rpad("n=$n", 12) for n in n_values]))
println("-" ^ 70)

for (name, ΔEL, R_load, R_drift, KT, S_trans, m0_trans) in results
    println(rpad(name, 38), join([rpad(round(R_load^n, digits=3), 12) for n in n_values]))
end

# ─────────────────────────────────────────────────────────────────────────────
# CYCLE COUNTING — RAYLEIGH AMPLITUDE BINS
# ─────────────────────────────────────────────────────────────────────────────
T_year   = 365.25 * 24 * 3600   # seconds per year
T_mean   = 2π / ω_p             # mean wave period, narrow-band approximation
N_cycles = T_year / T_mean      # wave cycles per year

n_bins = 200
a_max  = 4.0 * η0_p             # 4× peak wave amplitude, covers Rayleigh tail
a_bins = range(0, a_max, length=n_bins+1)
da     = (a_max - 0) / n_bins  # instead of step(a_bins)
a_mid  = collect(a_bins[1:end-1] .+ da/2)

println()
println("=" ^ 70)
println("CYCLE COUNT REDUCTION — amplitude bins, design life $T_design yr")
println("=" ^ 70)
println("(ratio is independent of annual exposure fraction)")
println()

for (name, ΔEL, R_load, R_drift, KT, S_trans, m0_trans) in results
    σ_inc   = sqrt(m0_inc)
    σ_trans = sqrt(m0_trans)

    p_inc   = (a_mid ./ σ_inc^2)   .* exp.(-a_mid.^2 ./ (2σ_inc^2))
    p_trans = (a_mid ./ σ_trans^2) .* exp.(-a_mid.^2 ./ (2σ_trans^2))

    N_inc   = N_cycles .* p_inc   .* da
    N_trans = N_cycles .* p_trans .* da

    tail_mask        = a_mid .> η0_p
    cycles_elim_tail = sum((N_inc .- N_trans)[tail_mask]) * T_design

    println("  $name")
    println("    σ_inc   = $(round(σ_inc,   digits=4)) m")
    println("    σ_trans = $(round(σ_trans, digits=4)) m")
    @printf("    Cycles above η0_p eliminated over %d yr: %.3e\n", T_design, cycles_elim_tail)
    println()
end