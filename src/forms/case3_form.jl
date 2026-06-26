# case3_form.jl
#
# Defines the frequency-domain weak form for Case 3: potential flow
# coupled to a rigid pipe's heave and pitch motion via hydrostatic and
# radiation coupling terms.

"Builds the Case 3 frequency-domain system: assembles the fluid weak
form (potential flow + damped free surface, no structure) together
with the fluid-rigid coupling blocks and hydrostatic stiffness for a
pipe free to heave and pitch, into one monolithic block matrix."
function build_case3_operator(sp, reg, meas, order, params)

  # -- Parameters -----------------------------------------------------------
  (; H, ρw, g, η₀, k, ω, μ₀, x₀, xdᵢₙ, xdₒᵤₜ, xc, zc, Mrb, Crb, Krb) = params

  # -- Incident wave fields -------------------------------------------------
  ηᵢₙ(x) = η₀ * exp(im*k*x[1])
  ϕᵢₙ(x) = -im*(η₀*ω/k) * (cosh(k*(x[2] + H)) / sinh(k*H)) * exp(im*k*x[1])

  uᵢₙ(x)  = (η₀*ω) * (cosh(k*(x[2] + H)) / sinh(k*H)) * exp(im*k*x[1])
  uzᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1])
  uₒᵤₜ(x) = 0.0 + 0im

  # -- Damping functions ----------------------------------------------------
  Ld = xdᵢₙ - x₀
  μ₁ᵢₙ(x) = μ₀ * (1.0 - sin(π/2 * (x[1]) / Ld))
  μ₁ₒᵤₜ(x) = μ₀ * (1.0 - cos(π/2 * (x[1] - xdₒᵤₜ) / Ld))
  μ₂ᵢₙ(x) = μ₁ᵢₙ(x) * k
  μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x) * k
  ηd(x)    = μ₂ᵢₙ(x) * ηᵢₙ(x)
  ∇ₙϕd(x)  = μ₁ᵢₙ(x) * uzᵢₙ(x)

  # -- Stabilization parameters ----------------------------------------------
  βh = 0.5
  αf = -im * ω / g

  # -- Measures --------------------------------------------------------------
  (; dΩ, dΓfs_mid, dΓd1, dΓd2, dΓpipe, dΓinlet, dΓoutlet) = meas
  (; npipe) = reg

  # -- Upward normal derivative used on free surface, same as Case 1 ---------
  nfs = reg.nfs
  ∇ₙfs(ϕ) = ∇(ϕ) ⋅ nfs

  # -- Rigid-body normal mode shapes on Γpipe --------------------------------
  ψ₁(x) = VectorValue(0.0, 1.0)                     # heave
  ψ₂(x) = VectorValue(-(x[2] - zc), x[1] - xc)      # pitch

  H1 = npipe ⋅ ψ₁   # heave normal displacement
  H2 = npipe ⋅ ψ₂   # pitch normal displacement

  # -- FE bilinear form for (ϕ, κ), Case 1 fluid/free-surface part, without beam terms --
  a_fe((ϕ, κ), (w, s)) = begin
    t1 = ∫( ∇(w) ⋅ ∇(ϕ) )dΩ
    t2 = ∫( βh * (s + αf*w) * (-im*ω*ϕ + g*κ) + im*ω*w*κ )dΓfs_mid
    t3 = ∫( βh * (s + αf*w) * (-im*ω*ϕ + g*κ) + im*ω*w*κ - μ₂ᵢₙ * κ * w + μ₁ᵢₙ * ∇ₙfs(ϕ) * (s + αf*w) )dΓd1
    t4 = ∫( βh * (s + αf*w) * (-im*ω*ϕ + g*κ) + im*ω*w*κ - μ₂ₒᵤₜ * κ * w + μ₁ₒᵤₜ * ∇ₙfs(ϕ) * (s + αf*w) )dΓd2
    t1 + t2 + t3 + t4
  end

  l_fe((w, s)) =
      ∫( w * uᵢₙ )dΓinlet +
      ∫( w * uₒᵤₜ )dΓoutlet -
      ∫( ηd * w - ∇ₙϕd * (s + αf*w) )dΓd1

  op_fe = AffineFEOperator(a_fe, l_fe, sp.X, sp.Y)
  Aff = get_matrix(op_fe)
  bfe = get_vector(op_fe)

  # -- Fluid <- rigid coupling block A_fξ: comes from iω ∫_Γpipe w ηp dΓ -------
  a_fξ1((w, s)) = ∫( im * ω * w * H1 )dΓpipe  # heave contribution
  a_fξ2((w, s)) = ∫( im * ω * w * H2 )dΓpipe  # pitch contribution

  Afξ1 = assemble_vector(a_fξ1, sp.Y)
  Afξ2 = assemble_vector(a_fξ2, sp.Y)
  Afξ  = hcat(Afξ1, Afξ2)     # puts the two coupling vectors together into a matrix with 2 columns (heave, pitch) and nfe rows

  # -- Rigid <- fluid coupling block A_ξf: comes from iωρw ∫_Γpipe δηp ϕ dΓ --
  a_ξf1((ϕ, κ)) = ∫( im * ω * ρw * H1 * ϕ )dΓpipe # heave contribution
  a_ξf2((ϕ, κ)) = ∫( im * ω * ρw * H2 * ϕ )dΓpipe # pitch contribution

  Aξf1 = assemble_vector(a_ξf1, sp.X)
  Aξf2 = assemble_vector(a_ξf2, sp.X)
  Aξf  = vcat(permutedims(Aξf1), permutedims(Aξf2)) # puts the two coupling vectors together into a matrix with 2 rows (heave, pitch) and nfe columns

  # -- Hydrostatic contribution on Γpipe: -ρw g ∫ δηp ηp dΓ ------------------
  # Hydrostatic stiffness matrix entries
  Khs11 = sum(∫( ρw * g * H1 * H1 )dΓpipe)
  Khs12 = sum(∫( ρw * g * H1 * H2 )dΓpipe)
  Khs21 = sum(∫( ρw * g * H2 * H1 )dΓpipe)
  Khs22 = sum(∫( ρw * g * H2 * H2 )dΓpipe)

  Khs = ComplexF64[
    Khs11 Khs12
    Khs21 Khs22
  ]

  # -- M ξ¨ + C ξ̇ + K ξ = Fhyd, frequency domain => (-ω² M - iω C + K) --------
  M = Mrb
  C = Crb
  K = Krb

  Arb = -ω^2 * M - im*ω * C + K - Khs

  # -- Combine into monolithic operator --------------------------------------
  A = [Aff  Afξ
       Aξf  Arb]

  b = vcat(bfe, zeros(ComplexF64, 2))

  return (A=A, b=b, Aff=Aff, bfe=bfe, Afξ=Afξ, Aξf=Aξf, Khs=Khs, Arb=Arb)
end