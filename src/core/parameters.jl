# parameters.jl
#
# Builds the physical and geometric parameter sets used by each case.
# Each function returns a NamedTuple of parameters, either extracted
# from the mesh's tagged boundary coordinates (Cases 1-3) or hardcoded
# to match a published reference geometry (Khabakhpasheva benchmark).

using Gridap.Geometry: get_node_coordinates

"Hardcoded parameter set matching the Khabakhpasheva et al. benchmark geometry, 
for validating Case 1 against the published reference case."
function build_khabakhpasheva_params()
  # Geometry 
  Lb = 12.5
  Ld = Lb
  Lf = 25.0
  LΩ = 2 * Ld + Lf          # = 50
  H  = 1.1               # water depth [m]
  α = 0.249              # Wavelength-to-beam ratio

  x₀ = 0.0

  xdᵢₙ = x₀ + Ld
  xb₀ = xdᵢₙ + Lb/2
  xb₁ = xb₀ + Lb
  xdₒᵤₜ = LΩ - Ld

  # Beam properties
  ρw = 1025.0            # water density [kg/m^3]
  ρb = 900.0             # beam density [kg/m^3]
  hb = 0.2               # beam thickness [m]
  g  = 9.81              # gravity [m/s^2]
  EI₁ = 47100           # beam bending stiffness [N m^2] for the first 20% of the beam
  EI₂ = 471              # beam bending stiffness [N m^2] for the remaining 80% of the beam
  EI_fun = x -> x[1] < (xb₀ + 0.2*Lb) ? EI₁ : EI₂ #471.0 * (1000.0 * (x[1]<(xb₀ + 0.2*Lb)) + 1.0*(x[1]>=(xb₀ + 0.2*Lb)))           # beam bending stiffness [N m^2]
  d₀ = 8.1561e-3         # Draft from Khabakhpasheva (for the test)
  m  = ρw * d₀           # beam mass per unit length [kg/m]

  # Wave properties
  η₀ = 0.01              # incident wave amplitude [m]
  λ  = α * Lb            # wavelength [m]
  k  = 2π / λ
  ω  = sqrt(g*k*tanh(k*H))
  T  = 2π / ω

  # Damping
  μ₀ = 2.5

  return (
    Lb=Lb, 
    Ld=Ld, 
    Lf=Lf, 
    LΩ=LΩ, 
    H=H, 
    x₀=x₀,
    xdᵢₙ=xdᵢₙ,
    xb₀=xb₀,
    xb₁=xb₁,
    xdₒᵤₜ=xdₒᵤₜ,
    α=α,
    ρw=ρw, 
    ρb=ρb, 
    hb=hb, 
    g=g, 
    EI₁=EI₁, 
    EI₂=EI₂,
    EI_fun=EI_fun, 
    d₀=d₀, 
    m=m,
    η₀=η₀, 
    λ=λ, 
    k=k, 
    ω=ω, 
    T=T,
    μ₀=μ₀
  )
end

"Builds the Case 1 parameter set: beam and wave properties for the thin elastic beam, 
with geometry extracted from the mesh's tagged boundaries."
function build_case1_params(model)
    # ── Extract geometry from mesh tags ──────────────────────────────────────
    grid     = get_grid(model)
    coords   = get_node_coordinates(grid)
    labeling = get_face_labeling(model)
    

    function coord_of_tag(tag)
        face_to_tag = labeling.d_to_dface_to_entity[1]
        tag_index   = findfirst(==(tag), labeling.tag_to_name)
        entity_id   = labeling.tag_to_entities[tag_index][1]
        node_ids    = findall(==(entity_id), face_to_tag)
        return coords[node_ids[1]]
    end

    xb₀   = coord_of_tag("BeamLeft")[1]
    xb₁   = coord_of_tag("BeamRight")[1]
    xdᵢₙ  = coord_of_tag("DampingInlet")[1]
    xdₒᵤₜ = coord_of_tag("DampingOutlet")[1]
    H     = abs(coord_of_tag("BottomLeft")[2])

    Lb  = xb₁ - xb₀
    Ld  = xdᵢₙ
    Lf  = xdₒᵤₜ - xdᵢₙ
    LΩ  = 2.0 * Ld + Lf
    x₀  = 0.0

    # ── Physical parameters ───────────────────────────────────────────────────
    ρw  = 1000.0
    ρb  = 62.5      # bare island
    hb  = 0.20
    g   = 9.81
    η₀  = 0.01
    μ₀  = 2.5

    # ── EI for SDR11 and SDR17 (Willemspolder D = 0.45m) ─────────────────────
    EI_SDR17 = 53930.0   # [Nm²/m]
    EI_SDR11 = 75435.0   # [Nm²/m]
    EI       = EI_SDR17   # default

    # ── Beam mass per unit length ─────────────────────────────────────────────
    m  = ρb * hb

    # ── Willemspolder wave parameters (peak frequency) ────────────────────────
    ωp = 2.581
    kp = 0.679
    λ  = 2π / kp
    k  = kp
    ω  = ωp
    T  = 2π / ω

    return (
        Lb=Lb, Ld=Ld, Lf=Lf, LΩ=LΩ, H=H, x₀=x₀,
        xb₀=xb₀, xb₁=xb₁, xdᵢₙ=xdᵢₙ, xdₒᵤₜ=xdₒᵤₜ,
        ρw=ρw, ρb=ρb, hb=hb, g=g, η₀=η₀, μ₀=μ₀,
        EI=EI, EI_SDR11=EI_SDR11, EI_SDR17=EI_SDR17,
        m=m,
        λ=λ, k=k, kp=kp, ω=ω, ωp=ωp, T=T
    )
end

"Builds the Case 2 parameter set for the time-domain formulation: beam, wave, and porous root-zone properties, 
with geometry extracted from the mesh's tagged boundaries. αₚ sets the porous resistance value for the sweep."
function build_case2_transient_params(model, αₚ = 0.0)
    grid     = get_grid(model)
    coords   = get_node_coordinates(grid)
    labeling = get_face_labeling(model)

    function coord_of_tag(tag)
        face_to_tag = labeling.d_to_dface_to_entity[1]
        tag_index   = findfirst(==(tag), labeling.tag_to_name)
        entity_id   = labeling.tag_to_entities[tag_index][1]
        node_ids    = findall(==(entity_id), face_to_tag)
        return coords[node_ids[1]]
    end

    xdᵢₙ  = coord_of_tag("DampingInlet")[1]
    xdₒᵤₜ = coord_of_tag("DampingOutlet")[1]
    xb₀   = coord_of_tag("BeamLeft")[1]
    xb₁   = coord_of_tag("BeamRight")[1]
    H = coord_of_tag("FreeSurfaceLeft")[2]   # z=0 at bottom, so H = free surface z

    x₀  = 0.0
    Ld  = xdᵢₙ
    Lf  = xdₒᵤₜ - xdᵢₙ
    LΩ  = 2*Ld + Lf
    Lb  = xb₁ - xb₀

    # ── Physical parameters ───────────────────────────────────────────────────
    ρw  = 1000.0
    ρb  = 62.5      # bare island
    hb  = 0.20
    g   = 9.81
    η₀  = 0.01
    ν   = 1e-3
    m   = ρb * hb

    # ── Beam geometry ─────────────────────────────────────────────────────────
    d₀  = ρb*hb/ρw     # draft [m]

    # ── EI ────────────────────────────────────────────────────────────────────
    EI_SDR17 = 53_930.0
    EI_SDR11 = 75_435.0
    EI       = EI_SDR17

    # ── Wave parameters ───────────────────────────────────────────────────────
    ωp = 2.581
    kp = 0.679
    λ  = 2π / kp
    k  = kp
    ω  = ωp
    T  = 2π / ω

    # ── Nitsche and porosity ──────────────────────────────────────────────────
    γ₀  = 10.0
    αₚ  = αₚ  # porosity drag value, sweep value


    return (
        Ld=Ld, Lf=Lf, LΩ=LΩ, H=H, x₀=x₀,
        xdᵢₙ=xdᵢₙ, xdₒᵤₜ=xdₒᵤₜ,
        xb₀=xb₀, xb₁=xb₁, Lb=Lb,
        ρw=ρw, ρb=ρb, hb=hb, g=g, η₀=η₀, ν=ν,
        EI=EI, EI_SDR11=EI_SDR11, EI_SDR17=EI_SDR17,
        m=m, d₀=d₀, 
        λ=λ, k=k, kp=kp, ω=ω, ωp=ωp, T=T,
        γ₀=γ₀, αₚ=αₚ
    )
end

"Builds the Case 3 parameter set: rigid pipe geometry derived from the mesh's pipe boundary, 
plus mass, inertia, and mooring matrices for the rigid-body heave/pitch response."
function build_case3_params(model)
  # ------------------------------------------------------------------
  # Extract geometry from mesh tags
  # ------------------------------------------------------------------
  grid     = get_grid(model)
  coords   = get_node_coordinates(grid)
  labeling = get_face_labeling(model)

  function coord_of_tag(tag)
    face_to_tag = labeling.d_to_dface_to_entity[1]
    tag_index   = findfirst(==(tag), labeling.tag_to_name)
    entity_id   = labeling.tag_to_entities[tag_index][1]
    node_ids    = findall(==(entity_id), face_to_tag)
    return coords[node_ids[1]]
  end

  xdᵢₙ  = coord_of_tag("DampingInlet")[1]
  xdₒᵤₜ = coord_of_tag("DampingOutlet")[1]
  H     = abs(coord_of_tag("BottomLeft")[2])

  # Find pipe arc nodes via Γpipe triangulation
  Γpipe_temp = BoundaryTriangulation(model, tags=["Pipe"])
  Γpipe_coords = get_cell_coordinates(Γpipe_temp)
  ncells = num_cells(Γpipe_temp)

  pipe_x = Float64[]
  pipe_z = Float64[]
  for i in 1:ncells
      for pt in Γpipe_coords[i]
          push!(pipe_x, pt[1])
          push!(pipe_z, pt[2])
      end
  end

  # Derive pipe geometry from bottom node
  z_bottom  = minimum(pipe_z)     # = zc - R = -draft
  draft     = abs(z_bottom)
  D         = 3.0 * draft / 2.0   # assumed draft-to-diameter ratio for a floating HDPE pipe (filled 2/3 with water)
  Rₒᵤₜ       = D / 2.0
  zc        = z_bottom + Rₒᵤₜ
  D         = 2.0 * Rₒᵤₜ
  Ld        = xdᵢₙ
  Lf        = xdₒᵤₜ - xdᵢₙ
  LΩ        = 2.0 * Ld + Lf
  xc        = LΩ / 2
  xpL       = minimum(pipe_x[pipe_z .≈ 0])  # waterline intersections at z≈0
  xpR       = maximum(pipe_x[pipe_z .≈ 0])
  xoff      = xpR - xc
  x₀        = 0.0

  # ------------------------------------------------------------------
  # Hardcoded physical parameters
  # ------------------------------------------------------------------
  t    = 0.01                 # pipe wall thickness [m]
  Rᵢₙ  = Rₒᵤₜ - t            # inner radius [m]
  ρp   = 930.0                # HDPE density [kg/m³]
  ρw   = 1000.0               # water density [kg/m³]
  ρair = 1.225                # air density [kg/m³]
  g    = 9.81                 # gravity [m/s²]
  η₀   = 0.01                 # incident wave amplitude [m]
  μ₀   = 2.5                  # damping coefficient [-]

  # ------------------------------------------------------------------
  # Willemspolder wave parameters (peak frequency)
  # ------------------------------------------------------------------
  ωp = 2.581                  # peak angular frequency [rad/s]
  kp = 0.679                  # peak wavenumber [rad/m]
  λ  = 2π / kp
  k  = kp
  ω  = ωp
  T  = 2π / ω

  # ------------------------------------------------------------------
  # Matrices inputs
  # ------------------------------------------------------------------
  m = ρp * π * (Rₒᵤₜ^2 - Rᵢₙ^2) + ( 2/3 * ρw * π * Rᵢₙ^2 ) + ( 1/3 * ρair * π * Rᵢₙ^2 )  # mass per unit length of the pipe [kg/m]
  Iwall = 0.5 * ρp * π * (Rₒᵤₜ^4 - Rᵢₙ^4)
  Iinside = 0.5 * ( (2/3) * ρw + (1/3) * ρair ) * π * Rᵢₙ^4
  Iₚ = Iwall + Iinside      # mass moment of inertia per unit length of the pipe [kg ⋅ m] 

  # Matrices for rigid-body motion
  Mrb = ComplexF64[
    m   0.0
    0.0 Iₚ
  ]

  # Mooring parameters
  k_heave = 0.0   # mooring stiffness in heave [N/m] — sweep over this
  k_pitch = 0.0   # mooring stiffness in pitch [Nm/rad] — sweep over this
  c_heave = 0.0   # mooring damping in heave [Ns/m] — sweep over this
  c_pitch = 0.0   # mooring damping in pitch [Nms/rad] — sweep over this

  Crb = ComplexF64[
      c_heave  0.0
      0.0      c_pitch
  ]

  Krb = ComplexF64[
      k_heave  0.0
      0.0      k_pitch
  ]

  return (
    Ld=Ld,
    Lf=Lf,
    LΩ=LΩ,
    H=H,
    x₀=x₀,
    xdᵢₙ=xdᵢₙ,
    xdₒᵤₜ=xdₒᵤₜ,
    D=D,
    Rᵢₙ=Rᵢₙ,
    Rₒᵤₜ=Rₒᵤₜ,
    draft=draft,
    xc=xc,
    zc=zc,
    xoff=xoff,
    xpL=xpL,
    xpR=xpR,
    ρw=ρw,
    ρp=ρp,
    g=g,
    η₀=η₀,
    λ=λ,
    k=k,
    kp=kp,
    ω=ω,
    ωp=ωp,
    T=T,
    μ₀=μ₀,
    Mrb=Mrb,
    Crb=Crb,
    Krb=Krb, 
    k_heave=k_heave,
    k_pitch=k_pitch,
    c_heave=c_heave,
    c_pitch=c_pitch
  )
end