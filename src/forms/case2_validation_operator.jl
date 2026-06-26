# case2_validation_operator.jl
#
# Defines the time-domain residual and Jacobians used to validate
# Case 2 against the Khabakhpasheva et al. benchmark.

"Builds the Case 2 time-domain operator for the Khabakhpasheva benchmark validation"
function build_case2_khab_transient_operator(sp, reg, meas, order, params)

    (; g, ρw, αₚ, γ₀, EI_fun, m) = params

    # -- Penalty ------------------------------------------------------------
    xΓstr = get_cell_coordinates(reg.Γstr)
    hsum  = 0.0
    for xs in xΓstr
        hsum += norm(xs[2] - xs[1])
    end
    h_e = hsum / length(xΓstr)
    γ   = 1.0 * order * (order - 1) / h_e
    γN  = γ₀ / h_e

    (; dΩ, dΩp, dΓfs, dΓstr, dΛstr) = meas
    (; nfs, nstr, nΛstr) = reg

    # -- Spatially varying EI as CellField ------------------------------------
    EI = CellField(EI_fun, reg.Γstr)

    # -- Residual -------------------------------------------------------------
    res(t, (u,p,κ,η,ζ), (w,q,s,v,χ)) =
        ∫( (∂t(u) ⋅ w) )dΩ                                                  +
        ∫( -(∇⋅w) * p )dΩ                                                   +
        ∫( q * (∇⋅u) )dΩ                                                    +
        ∫( αₚ * (w ⋅ u) )dΩp                                                +
        ∫( g * κ * (w ⋅ nfs) )dΓfs                                          +
        ∫( g * s * (∂t(κ) - (u ⋅ nfs)) )dΓfs                                +
        ∫( p * (w ⋅ nstr) )dΓstr                                            +
        ∫( -q * (u ⋅ nstr - ζ) )dΓstr                                       +
        ∫( γN * (u ⋅ nstr - ζ) * (w ⋅ nstr) )dΓstr                          +
        ∫( v * (m/ρw * ∂t(ζ) + g*η - p) )dΓstr                             +
        ∫( EI/ρw * Δ(v) * Δ(η) )dΓstr                                       -
        ∫( EI/ρw * (mean(Δ(η)) * jump(∇(v) ⋅ nΛstr)                         +
                 mean(Δ(v)) * jump(∇(η) ⋅ nΛstr)) )dΛstr                    +
        ∫( EI/ρw * γ * jump(∇(v) ⋅ nΛstr) * jump(∇(η) ⋅ nΛstr) )dΛstr      +
        ∫( χ * (∂t(η) - ζ) )dΓstr

    # -- Spatial Jacobian -----------------------------------------------------
    jac(t, (u,p,κ,η,ζ), (du,dp,dκ,dη,dζ), (w,q,s,v,χ)) =
        ∫( -(∇⋅w) * dp )dΩ                                                  +
        ∫( q * (∇⋅du) )dΩ                                                   +
        ∫( αₚ * (w ⋅ du) )dΩp                                               +
        ∫( g * dκ * (w ⋅ nfs) )dΓfs                                         +
        ∫( -g * s * (du ⋅ nfs) )dΓfs                                        +
        ∫( dp * (w ⋅ nstr) )dΓstr                                           +
        ∫( -q * (du ⋅ nstr - dζ) )dΓstr                                     +
        ∫( γN * (du ⋅ nstr - dζ) * (w ⋅ nstr) )dΓstr                        +
        ∫( v * (g*dη - dp) )dΓstr                                           +
        ∫( EI/ρw * Δ(v) * Δ(dη) )dΓstr                                      -
        ∫( EI/ρw * (mean(Δ(dη)) * jump(∇(v) ⋅ nΛstr) +
                 mean(Δ(v))  * jump(∇(dη) ⋅ nΛstr)) )dΛstr                  +
        ∫( EI/ρw * γ * jump(∇(v) ⋅ nΛstr) * jump(∇(dη) ⋅ nΛstr) )dΛstr     +
        ∫( -χ * dζ )dΓstr

    # -- Time Jacobian --------------------------------------------------------
    jac_t(t, (u,p,κ,η,ζ), (dut,dpt,dκt,dηt,dζt), (w,q,s,v,χ)) =
        ∫( (dut ⋅ w) )dΩ                                                    +
        ∫( g * s * dκt )dΓfs                                                 +
        ∫( v * m/ρw * dζt )dΓstr                                            +
        ∫( χ * dηt )dΓstr

    return TransientFEOperator(res, jac, jac_t, sp.X, sp.Y)
end