# case1_khab_validation_operator.jl
#
# Defines the frequency-domain weak form used to validate Case 1
# against the Khabakhpasheva et al. benchmark. 

"Builds the weak form operator for the Khabakhpasheva benchmark validation"
function build_khabakhpasheva_operator(sp, reg, meas, order, params)

  # Parameters
  (; Ld, H, ρw, g, m, η₀, k, ω, μ₀, xdₒᵤₜ) = params

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
  
  # Free-surface condition coefficients
  βh = 0.5
  αf = -im*ω/g

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
  EI = CellField(params.EI_fun, reg.Γstr) # spatially-varying bending stiffness, matching the Khabakhpasheva reference geometry

  a((ϕ, κ, η), (w, s, v)) = begin
    t1 = ∫( ∇(w) ⋅ ∇(ϕ) )dΩ
    t2 = ∫( im * ω * w * η )dΓstr
    t3 = ∫( v * (((-ω^2 * m) + ρw * g) * η - im * ω * ρw * ϕ) )dΓstr
    t4 = ∫( EI * Δ(v) * Δ(η) )dΓstr
    t5 = ∫( EI * (mean(Δ(η)) * jump(∇(v) ⋅ nΛstr) + mean(Δ(v)) * jump(∇(η) ⋅ nΛstr)) )dΛstr
    t6 = ∫( EI * γ * jump(∇(v) ⋅ nΛstr) * jump(∇(η) ⋅ nΛstr))dΛstr
    t7 = ∫( βh * (s + αf * w) * (-im * ω * ϕ + g * κ) + im * ω * w * κ )dΓfs_mid
    t8 = ∫( βh * (s + αf * w) * (-im * ω * ϕ + g * κ) + im * ω * w * κ - μ₂ᵢₙ * κ * w + μ₁ᵢₙ * ∇ₙ(ϕ) * (s + αf * w) )dΓd1
    t9 = ∫( βh * (s + αf * w) * (-im * ω * ϕ + g * κ) + im * ω * w * κ - μ₂ₒᵤₜ * κ *w + μ₁ₒᵤₜ *∇ₙ(ϕ) * (s + αf * w) )dΓd2

    t1 + t2 + t3 + t4 - t5 + t6 + t7 + t8 + t9
  end

  # Linear form
  l((w, s, v)) = ∫( w * uᵢₙ )dΓinlet + ∫( w * uₒᵤₜ )dΓoutlet - ∫( ηd * w - ∇ₙϕd * (s + αf * w) )dΓd1

  return AffineFEOperator(a, l, sp.X, sp.Y)
end