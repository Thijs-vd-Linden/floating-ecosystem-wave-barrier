# spaces.jl
#
# Builds the finite element trial/test spaces for each case, from the
# regions defined in domain.jl. Case 1 and Case 3 use complex-valued
# (ComplexF64) spaces for the frequency-domain potential-flow
# formulation. Case 2's time-domain formulation uses real-valued,
# time-dependent (Transient) spaces instead, with Dirichlet boundary
# conditions prescribed as functions of both position and time.

"Builds the Case 1 trial/test spaces: velocity potential on the fluid domain, 
free-surface elevation, and beam displacement, as a monolithic multi-field space."
function build_case1_spaces(reg; order::Int, T::Type=ComplexF64)
  # Reference elements (real basis; complex values via vector_type)
  reffe_Ω = ReferenceFE(lagrangian, Float64, order)
  reffe_Γfs = ReferenceFE(lagrangian, Float64, order)
  reffe_Γstr = ReferenceFE(lagrangian, Float64, order)

  # Fluid space on Ω: ϕ ∈ H1(Ω)
  V_Ω = TestFESpace(reg.Ω, reffe_Ω; conformity = :H1, vector_type = Vector{T})
  U_Ω = TrialFESpace(V_Ω) 

  # Surface/structure space on Γtop: η ∈ H1(Γtop)
  V_κ = TestFESpace(reg.Γfs, reffe_Γfs; conformity=:H1, vector_type = Vector{T})
  V_η = TestFESpace(reg.Γstr, reffe_Γstr; conformity=:H1, vector_type = Vector{T})
  U_κ = TrialFESpace(V_κ)
  U_η = TrialFESpace(V_η)

  # Monolithic multi-field spaces: (ϕ, η)d
  X = MultiFieldFESpace([U_Ω, U_κ, U_η])
  Y = MultiFieldFESpace([V_Ω, V_κ, V_η])

  return (
    reffe_Ω=reffe_Ω, reffe_Γfs=reffe_Γfs, reffe_Γstr=reffe_Γstr, 
    V_Ω = V_Ω, U_Ω = U_Ω,
    V_κ = V_κ, U_κ = U_κ, 
    V_η = V_η, U_η = U_η, 
    X = X, Y = Y)
end

"Builds the Case 2 time-domain trial/test spaces: velocity, pressure, free-surface elevation, beam displacement, and beam velocity, as a transient multi-field space. 
Inlet and free-surface Dirichlet conditions are prescribed as time-dependent incident wave functions."
function build_case2_transient_spaces(reg, params; order::Int)

  (; η₀, ω, k, H, ρw, ρb, hb, g) = params

  reffe_u    = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
  reffe_p    = ReferenceFE(lagrangian, Float64, order-1)
  reffe_Γfs  = ReferenceFE(lagrangian, Float64, order)
  reffe_Γstr = ReferenceFE(lagrangian, Float64, order)

  # ── Dirichlet functions (z∈[0,H] convention) ──────────────────────────────
  uin(x, t::Real) = VectorValue(
     η₀*ω * cosh(k*x[2]) / sinh(k*H) * cos(k*x[1] - ω*t),
    -η₀*ω * sinh(k*x[2]) / sinh(k*H) * sin(k*x[1] - ω*t)
  )
  uin(t::Real)    = x -> uin(x, t)

  u_bot(x, t::Real)  = VectorValue(0.0, 0.0)
  u_bot(t::Real)     = x -> u_bot(x, t)

  u_zero(x, t::Real) = VectorValue(0.0, 0.0)
  u_zero(t::Real)    = x -> u_zero(x, t)

  ηin(x, t::Real) = η₀ * cos(k*x[1] - ω*t)
  ηin(t::Real)    = x -> ηin(x, t)

  η_zero(x, t::Real) = 0.0
  η_zero(t::Real)    = x -> η_zero(x, t)

  # ── Velocity ──────────────────────────────────────────────────────────────
  V_u = TestFESpace(reg.Ω, reffe_u; conformity=:H1,
          dirichlet_tags=[TAG_inlet, TAG_outlet, TAG_bot],
          dirichlet_masks=[(true,true), (true,true), (true,true)])
  U_u = TransientTrialFESpace(V_u, [uin, u_zero, u_bot])

  # ── Pressure (C0, not transient) ──────────────────────────────────────────
  V_p = TestFESpace(reg.Ω, reffe_p; conformity=:L2)
  U_p = TrialFESpace(V_p)

  # ── Free surface elevation ─────────────────────────────────────────────────
  V_κ = TestFESpace(reg.Γfs, reffe_Γfs; conformity=:H1,
          dirichlet_tags=[TAG_fsLeft, TAG_fsRight])
  U_κ = TransientTrialFESpace(V_κ, [ηin, η_zero])

  # ── Beam displacement (free ends) ─────────────────────────────────────────
  V_η = TestFESpace(reg.Γstr, reffe_Γstr; conformity=:H1)
  U_η = TransientTrialFESpace(V_η)

  # ── Beam velocity ζ = ∂t(η) (free ends) ──────────────────────────────────
  V_ζ = TestFESpace(reg.Γstr, reffe_Γstr; conformity=:H1)
  U_ζ = TransientTrialFESpace(V_ζ)

  X = TransientMultiFieldFESpace([U_u, U_p, U_κ, U_η, U_ζ])
  Y = MultiFieldFESpace([V_u, V_p, V_κ, V_η, V_ζ])

  return (
    V_u=V_u, U_u=U_u,
    V_p=V_p, U_p=U_p,
    V_κ=V_κ, U_κ=U_κ,
    V_η=V_η, U_η=U_η,
    V_ζ=V_ζ, U_ζ=U_ζ,
    X=X, Y=Y
  )
end

"Builds the Case 3 trial/test spaces: velocity potential on the fluid domain and free-surface elevation, as a monolithic multi-field space."
function build_case3_spaces(reg, params; order::Int, T::Type=ComplexF64)
  # Reference elements (real basis; complex values via vector_type)
  reffe_Ω = ReferenceFE(lagrangian, Float64, order)
  reffe_Γfs = ReferenceFE(lagrangian, Float64, order)

  # Fluid space on Ω: ϕ ∈ H1(Ω)
  V_Ω = TestFESpace(reg.Ω, reffe_Ω; conformity = :H1, vector_type = Vector{T})
  U_Ω = TrialFESpace(V_Ω) 

  # Free-surface space on Γfs: κ ∈ H1(Γfs)
  V_κ = TestFESpace(reg.Γfs, reffe_Γfs; conformity=:H1, vector_type=Vector{T})
  U_κ = TrialFESpace(V_κ)

  # Monolithic FE spaces: (ϕ, κ)
  X = MultiFieldFESpace([U_Ω, U_κ])
  Y = MultiFieldFESpace([V_Ω, V_κ])

  return (
    reffe_Ω=reffe_Ω, reffe_Γfs=reffe_Γfs,
    V_Ω=V_Ω, U_Ω=U_Ω, 
    V_κ=V_κ, U_κ=U_κ,
    X=X, Y=Y)
end