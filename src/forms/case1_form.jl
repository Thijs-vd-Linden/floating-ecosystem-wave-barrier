# case1_form.jl
#
# Defines the frequency-domain weak form for Case 1: potential flow
# coupled to a thin elastic beam via the free surface. The structure
# is unmoored; mooring is only modelled for Case 3.

"Builds the Case 1 frequency-domain operator: the monolithic weak form coupling fluid potential flow, 
free-surface elevation, and the thin elastic beam (Euler-Bernoulli, free ends), with Sommerfeld-type damping zones at the inlet and outlet boundaries."
function build_case1_operator(sp, reg, meas, order, params)
  # Parameters
  (; Ld, H, ρw, ρb, hb, g, EI, η₀, k, ω, T, μ₀, xdₒᵤₜ) = params

  # Incident wave fields
  ηᵢₙ(x) = η₀ * exp(im*k*x[1])
  ϕᵢₙ(x) = -im*(η₀*ω/k) * (cosh(k*(x[2] + H)) / sinh(k*H)) * exp(im*k*x[1])

  uᵢₙ(x) = (η₀*ω) * (cosh(k*(x[2] + H)) / sinh(k*H)) * exp(im*k*x[1])
  uzᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1])
  uₒᵤₜ(x) = 0.0 + 0im

  # Damping functions
  μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1])/Ld))
  μ₁ₒᵤₜ(x) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒᵤₜ)/Ld))
  μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
  μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
  ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
  ∇ₙϕd(x) = μ₁ᵢₙ(x)*uzᵢₙ(x)

  # Stabilization parameters
  βh = 0.5
  αf = -im * ω / g

  # Penalty parameter: now based on Γstr, not Γtop
  xΓstr = get_cell_coordinates(reg.Γstr)
  hsum = 0.0
  for xs in xΓstr
    hsum += norm(xs[2] - xs[1])
  end
  h = hsum / length(xΓstr)
  γ = 1.0 * order * (order - 1) / h

    # Measures
  (; dΩ, dΓfs_mid, dΓd1, dΓd2, dΓstr, dΓinlet, dΓoutlet, dΛstr) = meas

  # Skeleton normal on structure only
  nΛstr = reg.nΛstr

  # Weak form
   ∇ₙ(ϕ) = ∇(ϕ) ⋅ VectorValue(0.0, 1.0)  # Normal derivative ∂ϕ/∂n 

  a((ϕ, κ, η), (w, s, v)) = begin
    t1 = ∫( ∇(w) ⋅ ∇(ϕ) )dΩ                                              # Laplace (fluid interior)
    t2 = ∫( im * ω * w * η )dΓstr                                        # kinematic FS-structure coupling
    t3 = ∫( v * (((-ω^2 * ρb * hb) + ρw * g) * η - im * ω * ρw * ϕ) )dΓstr  # beam inertia + hydrostatic restoring, driven by fluid pressure
    t4 = ∫( EI * Δ(v) * Δ(η) )dΓstr                                       # beam bending 
    t5 = ∫( EI * (mean(Δ(η)) * jump(∇(v) ⋅ nΛstr) + mean(Δ(v)) * jump(∇(η) ⋅ nΛstr)) )dΛstr   # beam bending (interface consistency)
    t6 = ∫( EI * γ * jump(∇(v) ⋅ nΛstr) * jump(∇(η) ⋅ nΛstr))dΛstr          # beam bending (penalty term)
    t7 = ∫( βh * (s + αf * w) * (-im * ω * ϕ + g * κ) + im * ω * w * κ )dΓfs_mid  # kinematic + dynamic FS condition
    t8 = ∫( βh * (s + αf * w) * (-im * ω * ϕ + g * κ) + im * ω * w * κ - μ₂ᵢₙ * κ * w + μ₁ᵢₙ * ∇ₙ(ϕ) * (s + αf * w) )dΓd1   # FS condition + inlet damping zone
    t9 = ∫( βh * (s + αf * w) * (-im * ω * ϕ + g * κ) + im * ω * w * κ - μ₂ₒᵤₜ * κ *w + μ₁ₒᵤₜ *∇ₙ(ϕ) * (s + αf * w) )dΓd2    # FS condition + outlet damping zone

    t1 +t2 + t3 + t4 - t5 + t6 + t7 + t8 + t9
  end

  # Linear form: inlet forcing | outlet forcing | incident-field correction for the inlet damping zone
  l((w, s, v)) = ∫( w * uᵢₙ )dΓinlet + ∫( w * uₒᵤₜ )dΓoutlet - ∫( ηd * w - ∇ₙϕd * (s + αf * w) )dΓd1

  return AffineFEOperator(a, l, sp.X, sp.Y)
end