# case2_sweep_L18.jl
#
# Time-domain porous-resistance sweep for Case 2, beam length L_b =
# 18.4 m. Identical structure to case2_sweep_L4.jl, but uses the
# longer L_18 domain rather than a damping zone to limit outlet-
# boundary reflection from contaminating the extraction window.
#
# Run case2_time_postprocess.jl afterwards to reanalyse the saved
# JLD2 output and produce the figures used in the thesis.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include("../src/Thesis_FPVBarrier.jl")
using .Thesis_FPVBarrier
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Visualization: paraview_collection, vtk_save, createvtk
using Plots
using CSV
using DataFrames
using JLD2


# Plot function
function make_plots(ts, κ_profiles, κ_max_t, η_max_t,
                    κ_T1_t, κ_R1_t,
                    xs_κ, xs_κ_raw, i_T1, i_R1,
                    x_T1, x_T2, x_T3, x_R1, x_R2, x_R3,
                    T, LΩ, xdᵢₙ, xdₒᵤₜ, xb₀, xb₁, η₀,
                    win_start, win_end, inc_start1, inc_end1, inc_end2,
                    αp, ap_str, plot_dir, prefix)

    # Max amplitude over time
    plt_max = plot(ts ./ T, κ_max_t,
        label="max |κ| / η₀", xlabel="t / T", ylabel="max |κ| / η₀",
        title="Max free surface amplitude — αp = $(αp) s⁻¹",
        lw=1.5, framestyle=:box, grid=true)
    vline!(plt_max, [inc_start1], ls=:dot, color=:purple, label="inc window start")
    vline!(plt_max, [inc_end2],   ls=:dot, color=:purple, label="inc window end")
    vline!(plt_max, [win_start],  ls=:dot, color=:black,  label="extraction start")
    vline!(plt_max, [win_end],    ls=:dot, color=:gray,   label="extraction end")
    savefig(plt_max, joinpath(plot_dir, "$(prefix)_maxamp_ap$(ap_str).png"))

    # Beam max deflection over time
    plt_eta = plot(ts ./ T, η_max_t,
        label="max |η| / η₀", xlabel="t / T", ylabel="max |η| / η₀",
        title="Max beam deflection — αp = $(αp) s⁻¹",
        lw=1.5, framestyle=:box, grid=true)
    vline!(plt_eta, [win_start], ls=:dot, color=:black, label="extraction start")
    vline!(plt_eta, [win_end],   ls=:dot, color=:gray,  label="extraction end")
    savefig(plt_eta, joinpath(plot_dir, "$(prefix)_beam_ap$(ap_str).png"))

    # Free surface envelope over extraction window
    i_env = findall(tᵢ -> win_start*T <= tᵢ <= win_end*T, ts)
    κ_env = vec(maximum(abs.(hcat(κ_profiles[i_env]...)), dims=2))
    plt_env = plot(xs_κ ./ LΩ, κ_env,
        label="envelope [$(win_start)T, $(win_end)T]",
        xlabel="x / LΩ", ylabel="max |κ| / η₀",
        title="Free surface envelope — αp = $(αp) s⁻¹",
        lw=1.5, framestyle=:box, grid=true)
    hline!(plt_env, [1.0],  ls=:dash, color=:black,  label="η₀")
    vline!(plt_env, [xdᵢₙ  / LΩ], ls=:dash, color=:red,   label="probe region left")
    vline!(plt_env, [xdₒᵤₜ / LΩ], ls=:dash, color=:green, label="probe region right")
    vline!(plt_env, [xb₀   / LΩ], ls=:dot,  color=:blue,  label="beam left")
    vline!(plt_env, [xb₁   / LΩ], ls=:dot,  color=:blue,  label="beam right")
    savefig(plt_env, joinpath(plot_dir, "$(prefix)_envelope_ap$(ap_str).png"))

    # Snapshot profiles at 13T, 14T, 15T
    i_13 = argmin(abs.(ts .- 13.0*T))
    i_14 = argmin(abs.(ts .- 14.0*T))
    i_15 = argmin(abs.(ts .- 15.0*T))
    plt_fs = plot(xlabel="x / LΩ", ylabel="κ / η₀",
        title="Free surface profiles — αp = $(αp) s⁻¹",
        framestyle=:box, grid=true)
    plot!(plt_fs, xs_κ ./ LΩ, κ_profiles[i_13], label="t = 13T", lw=1.5)
    plot!(plt_fs, xs_κ ./ LΩ, κ_profiles[i_14], label="t = 14T", lw=1.5, ls=:dash)
    plot!(plt_fs, xs_κ ./ LΩ, κ_profiles[i_15], label="t = 15T", lw=1.5, ls=:dot)
    vline!(plt_fs, [xdᵢₙ  / LΩ], ls=:dash, color=:red,   label="probe region left")
    vline!(plt_fs, [xdₒᵤₜ / LΩ], ls=:dash, color=:green, label="probe region right")
    vline!(plt_fs, [xb₀   / LΩ], ls=:dot,  color=:blue,  label="beam left")
    vline!(plt_fs, [xb₁   / LΩ], ls=:dot,  color=:blue,  label="beam right")
    savefig(plt_fs, joinpath(plot_dir, "$(prefix)_profiles_ap$(ap_str).png"))

    # Probe time series
    plt_ts = plot(ts ./ T, κ_T1_t ./ η₀,
        label="κ_T1 (x=$(round(xs_κ_raw[i_T1],digits=1)) m)",
        xlabel="t / T", ylabel="κ / η₀",
        title="Probe time series — αp = $(αp) s⁻¹",
        lw=1.5, framestyle=:box, grid=true)
    plot!(plt_ts, ts ./ T, κ_R1_t ./ η₀,
        label="κ_R1 (x=$(round(xs_κ_raw[i_R1],digits=1)) m)",
        lw=1.5, ls=:dash)
    vline!(plt_ts, [inc_start1], ls=:dot, color=:purple, label="inc window start")
    vline!(plt_ts, [inc_end2],   ls=:dot, color=:purple, label="inc window end")
    vline!(plt_ts, [win_start],  ls=:dot, color=:black,  label="extraction start")
    vline!(plt_ts, [win_end],    ls=:dot, color=:gray,   label="extraction end")
    savefig(plt_ts, joinpath(plot_dir, "$(prefix)_timeseries_ap$(ap_str).png"))
end


# -- Warmup ------------------------------------------------------------------------------
println("Warming up...")
let
    order = 2
    model_wm  = Thesis_FPVBarrier.build_model("meshes/case2_length_sweep/case2_L_18.msh")
    params_wm = Thesis_FPVBarrier.build_case2_transient_params(model_wm, 0.0)
    reg_wm    = Thesis_FPVBarrier.build_case2_transient_regions(model_wm)
    meas_wm   = Thesis_FPVBarrier.build_case2_transient_measures(reg_wm; order=order)
    sp_wm     = Thesis_FPVBarrier.build_case2_transient_spaces(reg_wm, params_wm; order=order)
    op_wm     = Thesis_FPVBarrier.build_case2_transient_operator(sp_wm, reg_wm, meas_wm, order, params_wm)
    p0_wm     = interpolate_everywhere(x -> params_wm.g*(params_wm.H - x[2]), sp_wm.U_p)
    x_wm      = interpolate_everywhere([VectorValue(0.0,0.0), p0_wm, 0.0, 0.0, 0.0], sp_wm.X(0.0))
    ls_wm     = LUSolver()
    xht_wm    = solve(ThetaMethod(ls_wm, params_wm.T/2, 0.5), op_wm, x_wm, 0.0, params_wm.T/2)
    for _ in xht_wm end
end
println("Warmup complete.")


function run_case2_sweep_L18()
    # -- Settings -------------------------------------------------------------------------
    order      = 2
    dt_per_T   = 40        # timesteps per period
    tend_T     = 17        # total run duration [T] — clean window [12T,15T] + buffer
    win_start  = 12      # extraction window start [T]
    win_end    = 15      # extraction window end [T]
    inc_start1 = 5.0       # incident window 1 start [T]
    inc_end1   = 8.0       # incident window 1 end [T]
    inc_start2 = 5.0       # incident window 2 start [T]
    inc_end2   = 11.0      # incident window 2 end [T]
    vtk_every  = 10        # save every nth step in VTK window
    prefix     = "case2_L18"

    αp_sweep = [0.00, 0.01, 0.1, 0.5, 1.0, 2.5, 5.0]
    
    plot_dir    = "results/plots/case2_timedomain/sweep/L18"
    vtk_base    = "results/vtk/case2_sweep_L18"
    jld2_dir    = "results/jld2/case2/L18"
    mkpath(plot_dir)
    mkpath(jld2_dir)

    # -- Sweep ---------------------------------------------------------------------------
    for αp in αp_sweep
        ap_str = replace(string(round(αp, digits=2)), "." => "p")
        println("\n========================================")
        println("Running αp = ", αp, " s⁻¹  (Lb = 18.4 m)")
        println("========================================")

        # Output directories
        vtk_dir_run = joinpath(vtk_base, "ap_$(ap_str)")
        rm(vtk_dir_run, recursive=true, force=true)
        mkpath(joinpath(vtk_dir_run, "fluid"))
        mkpath(joinpath(vtk_dir_run, "fs"))
        mkpath(joinpath(vtk_dir_run, "str"))

        fname_fluid = joinpath(vtk_dir_run, "fluid", "fluid_alphap$(ap_str)")
        fname_fs    = joinpath(vtk_dir_run, "fs",    "fs_alphap$(ap_str)")
        fname_str   = joinpath(vtk_dir_run, "str",   "str_alphap$(ap_str)")

        # Build model and params
        model2  = Thesis_FPVBarrier.build_model("meshes/case2_length_sweep/case2_L_18.msh")
        params2 = Thesis_FPVBarrier.build_case2_transient_params(model2, αp)
        reg2    = Thesis_FPVBarrier.build_case2_transient_regions(model2)
        meas2   = Thesis_FPVBarrier.build_case2_transient_measures(reg2; order=order)
        sp2     = Thesis_FPVBarrier.build_case2_transient_spaces(reg2, params2; order=order)
        op2     = Thesis_FPVBarrier.build_case2_transient_operator(sp2, reg2, meas2, order, params2)

        (; T, η₀, LΩ, xdᵢₙ, xdₒᵤₜ, xb₀, xb₁, k, ω, g, H) = params2

        dt         = T / dt_per_T
        tend       = tend_T * T
        ls         = LUSolver()
        ode_solver = ThetaMethod(ls, dt, 0.5)

        # Probe setup — close to beam
        λ    = 2π / k

        x_T1 = xb₁ + 1.0*λ
        Δx_T = min(λ/4, 2.0*λ / 2.5)
        x_T2 = x_T1 + Δx_T
        x_T3 = x_T1 + 2*Δx_T

        x_R1 = xb₀ - 2.0*λ
        Δx_R = min(λ/4, 2.0*λ / 2.5)
        x_R2 = x_R1 - Δx_R
        x_R3 = x_R1 - 2*Δx_R

        xs_κ_raw = get_free_dof_values(interpolate_everywhere(x -> x[1], sp2.V_κ))
        sort_idx = sortperm(xs_κ_raw)
        xs_κ     = xs_κ_raw[sort_idx]

        x_η_all    = get_free_dof_values(interpolate_everywhere(x -> x[1], sp2.V_η))
        sort_idx_η = sortperm(x_η_all)
        xs_η       = x_η_all[sort_idx_η]

        i_T1 = argmin(abs.(xs_κ_raw .- x_T1));  i_T2 = argmin(abs.(xs_κ_raw .- x_T2))
        i_T3 = argmin(abs.(xs_κ_raw .- x_T3));  i_R1 = argmin(abs.(xs_κ_raw .- x_R1))
        i_R2 = argmin(abs.(xs_κ_raw .- x_R2));  i_R3 = argmin(abs.(xs_κ_raw .- x_R3))

        println("Probes T: x = ", round.([xs_κ_raw[i_T1], xs_κ_raw[i_T2], xs_κ_raw[i_T3]], digits=3))
        println("Probes R: x = ", round.([xs_κ_raw[i_R1], xs_κ_raw[i_R2], xs_κ_raw[i_R3]], digits=3))
        println("xdᵢₙ = ", xdᵢₙ, "  xb₀ = ", xb₀, "  xb₁ = ", xb₁, "  xdₒᵤₜ = ", xdₒᵤₜ)

        # Storage
        κ_T1_t = Float64[];  κ_T2_t = Float64[];  κ_T3_t = Float64[]
        κ_R1_t = Float64[];  κ_R2_t = Float64[];  κ_R3_t = Float64[]
        κ_max_t    = Float64[]
        η_max_t    = Float64[]
        κ_profiles = Vector{Vector{Float64}}()
        η_profiles = Vector{Vector{Float64}}()

        # Initial conditions
        p0_fun = interpolate_everywhere(x -> g*(H - x[2]), sp2.U_p)
        x₀     = interpolate_everywhere([VectorValue(0.0,0.0), p0_fun, 0.0, 0.0, 0.0], sp2.X(0.0))

        # VTK collections
        pvd_Ω    = paraview_collection(fname_fluid, append=false)
        pvd_Γfs  = paraview_collection(fname_fs,    append=false)
        pvd_Γstr = paraview_collection(fname_str,   append=false)

        # Single-stage solver
        step = 0
        steps_per_period = round(Int, T/dt)
        xht  = solve(ode_solver, op2, x₀, 0.0, tend)

        println("Starting time integration: tend = $(tend_T) T, dt = $(round(dt, digits=4)) s")
        @time for ((uh, ph, κh, ηh, ζh), t) in xht
            step += 1
            κ_dofs = get_free_dof_values(κh)
            η_dofs = get_free_dof_values(ηh)

            if mod(step, steps_per_period) == 0
                n_period = div(step, steps_per_period)
                κ_max    = maximum(abs.(κ_dofs)) / η₀
                η_max    = maximum(abs.(η_dofs)) / η₀
                println("t = ", n_period, " T  |  max |κ|/η₀ = ", round(κ_max, digits=3),
                        "  |  max |η|/η₀ = ", round(η_max, digits=3))
                if isnan(κ_max) || isnan(η_max)
                    println("NaN detected, stopping.")
                    break
                end
            end

            push!(κ_T1_t, κ_dofs[i_T1]);  push!(κ_T2_t, κ_dofs[i_T2]);  push!(κ_T3_t, κ_dofs[i_T3])
            push!(κ_R1_t, κ_dofs[i_R1]);  push!(κ_R2_t, κ_dofs[i_R2]);  push!(κ_R3_t, κ_dofs[i_R3])
            push!(κ_max_t, maximum(abs.(κ_dofs)) / η₀)
            push!(η_max_t, maximum(abs.(η_dofs)) / η₀)
            push!(κ_profiles, copy(κ_dofs[sort_idx] ./ η₀))
            push!(η_profiles, copy(η_dofs[sort_idx_η] ./ η₀))
            
            if win_start*T <= t <= win_end*T && mod(step, vtk_every) == 0
                pvd_Ω[t]    = createvtk(reg2.Ω,    fname_fluid * "_$(step)", cellfields=["u" => uh, "p" => ph])
                pvd_Γfs[t]  = createvtk(reg2.Γfs,  fname_fs    * "_$(step)", cellfields=["kappa" => κh])
                pvd_Γstr[t] = createvtk(reg2.Γstr, fname_str   * "_$(step)", cellfields=["eta" => ηh, "zeta" => ζh])
            end
        end

        ts = collect(range(dt, tend, length=length(κ_T1_t)))
        println("Time integration complete at t = ", round(ts[end]/T, digits=2), " T")

        vtk_save(pvd_Ω);  vtk_save(pvd_Γfs);  vtk_save(pvd_Γstr)
        println("VTK saved for αp = ", αp)

        # -- Save JLD2 --------------------------------------------------------------------
        jld2_path = joinpath(jld2_dir, "$(prefix)_ap$(ap_str).jld2")
        @save jld2_path κ_profiles xs_κ xs_κ_raw ts αp LΩ xb₀ xb₁ xdᵢₙ xdₒᵤₜ η₀ T ω k κ_T1_t κ_T2_t κ_T3_t κ_R1_t κ_R2_t κ_R3_t κ_max_t η_max_t x_T1 x_T2 x_T3 x_R1 x_R2 x_R3 i_T1 i_T2 i_T3 i_R1 i_R2 i_R3 η_profiles xs_η
        println("JLD2 saved: ", jld2_path)

        # -- Plots ------------------------------------------------------------------------
        make_plots(ts, κ_profiles, κ_max_t, η_max_t,
                   κ_T1_t, κ_R1_t,
                   xs_κ, xs_κ_raw, i_T1, i_R1,
                   x_T1, x_T2, x_T3, x_R1, x_R2, x_R3,
                   T, LΩ, xdᵢₙ, xdₒᵤₜ, xb₀, xb₁, η₀,
                   win_start, win_end, inc_start1, inc_end1, inc_end2,
                   αp, ap_str, plot_dir, prefix)
        println("Plots saved for αp = ", αp)
    end
end

@time run_case2_sweep_L18()
