# case2_transient_form.jl
#
# Defines the time-domain residual and Jacobians for Case 2: a
# velocity-pressure formulation of potential flow, coupled to a thin
# elastic beam with a submerged porous root zone (Darcy drag). Used
# with a transient (time-stepping) solver; no free-surface damping
# zones are used here, unlike the frequency-domain Cases 1 and 3.

"Builds the Case 2 time-domain operator: residual, spatial Jacobian,
and time Jacobian for the coupled fluid-beam-porous-root-zone system."
function build_case2_transient_operator(sp, reg, meas, order, params)

  (; g, ρw, ρb, hb, EI, αₚ, γ₀) = params

  # -- Penalty ----------------------------------------------------------
  xΓstr = get_cell_coordinates(reg.Γstr)
  hsum  = 0.0
  for xs in xΓstr
    hsum += norm(xs[2] - xs[1])
  end
  h_e = hsum / length(xΓstr)
  γ   = 1.0 * order * (order - 1) / h_e
  γN  =  γ₀ / h_e 

  # -- Measures and normals ------------------------------------------------
  (; dΩ, dΩp, dΓfs, dΓstr, dΛstr) = meas
  (; nfs, nstr, nΛstr) = reg

  # -- Residual -------------------------------------------------------------
  res(t, (u,p,κ,η,ζ), (w,q,s,v,χ)) =
    # -- Fluid momentum -----------------------------------------------------
    ∫( (∂t(u) ⋅ w) )dΩ                                                      +
    ∫( -(∇⋅w) * p )dΩ                                                       +
    ∫( q * (∇⋅u) )dΩ                                                        +
    # -- Porous drag --------------------------------------------------------
    ∫( αₚ * (w ⋅ u) )dΩp                                                    +
    # -- Free surface -------------------------------------------------------
    ∫( g * κ * (w ⋅ nfs) )dΓfs                                              +
    ∫( g * s * (∂t(κ) - (u ⋅ nfs)) )dΓfs                                    +
    # -- Nitsche FSI --------------------------------------------------------
    ∫( p * (w ⋅ nstr) )dΓstr                                               +
    ∫( -q * (u ⋅ nstr - ζ) )dΓstr                                          +
    ∫( γN * (u ⋅ nstr - ζ) * (w ⋅ nstr) )dΓstr                              +
    # -- Beam equation ------------------------------------------------------
    ∫( v * (ρb*hb/ρw *∂t(ζ) + g*η - p) )dΓstr                              +
    ∫( EI/ρw * Δ(v) * Δ(η) )dΓstr                                          -
    ∫( EI/ρw * (mean(Δ(η)) * jump(∇(v) ⋅ nΛstr)                             +
             mean(Δ(v)) * jump(∇(η) ⋅ nΛstr)) )dΛstr                        +
    ∫( EI/ρw * γ * jump(∇(v) ⋅ nΛstr) * jump(∇(η) ⋅ nΛstr) )dΛstr           +
    # -- Auxiliary: ζ = ∂t(η) -----------------------------------------------
    ∫( χ * (∂t(η) - ζ) )dΓstr

  # -- Spatial Jacobian -----------------------------------------------------
  jac(t, (u,p,κ,η,ζ), (du,dp,dκ,dη,dζ), (w,q,s,v,χ))                        =
    # -- Fluid momentum -----------------------------------------------------
    ∫( -(∇⋅w) * dp )dΩ                                                      +
    ∫( q * (∇⋅du) )dΩ                                                       +
    # -- Porous drag --------------------------------------------------------
    ∫( αₚ * (w ⋅ du) )dΩp                                                    +
    # -- Free surface -------------------------------------------------------
    ∫( g * dκ * (w ⋅ nfs) )dΓfs                                             +
    ∫( -g * s * (du ⋅ nfs) )dΓfs                                            +
    # -- Nitsche FSI -------------------------------------------------------- 
    ∫( dp * (w ⋅ nstr) )dΓstr                                               +
    ∫( -q * (du ⋅ nstr - dζ) )dΓstr                                         +
    ∫( γN * (du ⋅ nstr - dζ) * (w ⋅ nstr) )dΓstr                            +
    # -- Beam equation ------------------------------------------------------
    ∫( v * (g*dη - dp) )dΓstr                                               +
    ∫( EI/ρw * Δ(v) * Δ(dη) )dΓstr                                          -
    ∫( EI/ρw * (mean(Δ(dη)) * jump(∇(v) ⋅ nΛstr) +
             mean(Δ(v))  * jump(∇(dη) ⋅ nΛstr)) )dΛstr                       +
    ∫( EI/ρw * γ * jump(∇(v) ⋅ nΛstr) * jump(∇(dη) ⋅ nΛstr) )dΛstr           +
    # -- Auxiliary ---------------------------------------------------------- 
    ∫( -χ * dζ )dΓstr

  # -- Time Jacobian --------------------------------------------------------
  jac_t(t, (u,p,κ,η,ζ), (dut,dpt,dκt,dηt,dζt), (w,q,s,v,χ))                  =
    ∫( (dut ⋅ w) )dΩ                                                         +
    ∫( g * s * dκt )dΓfs                                                     +
    ∫( v * ρb*hb/ρw * dζt )dΓstr                                             +
    ∫( χ * dηt )dΓstr

  return TransientFEOperator(res, jac, jac_t, sp.X, sp.Y)
end