# khabakhpasheva_domain.jl
#
# Builds the mesh model, tagged regions, and integration measures for
# the Khabakhpasheva benchmark validation domain. Mirrors the structure
# of domain.jl.

using GridapGmsh

"Load a Gmsh .msh file into a Gridap discrete model"
function build_khabakhpasheva_model(filename)
    model = GmshDiscreteModel(filename)
    return model
end

"Building tagged regions for the khabakhpasheva benchmark test"    # Defines where cells and faces are located in the domain, and how they are connected
function build_regions_khabakhpasheva(model, params)
    Ω = Triangulation(model)                                
    Γbot    = BoundaryTriangulation(model, tags=TAG_bot)
    Γinlet  = BoundaryTriangulation(model, tags=TAG_inlet)
    Γoutlet = BoundaryTriangulation(model, tags=TAG_outlet)

    Γfs  = BoundaryTriangulation(model, tags=TAG_fs)
    Γstr = BoundaryTriangulation(model, tags=TAG_str)
    Λstr  = Skeleton(Γstr)
    nΛstr = get_normal_vector(Λstr)

    xΓfs = get_cell_coordinates(Γfs)

    function is_damping1(xs)
        n = length(xs)
        xmid = sum(x[1] for x in xs) / n
        return params.x₀ <= xmid <= params.xdᵢₙ
    end

    function is_damping2(xs)
        n = length(xs)
        xmid = sum(x[1] for x in xs) / n
        return params.xdₒᵤₜ <= xmid <= params.LΩ
    end

    Γd1_to_Γfs_mask = lazy_map(is_damping1, xΓfs)
    Γd2_to_Γfs_mask = lazy_map(is_damping2, xΓfs)

    Γd1_to_Γfs = findall(Γd1_to_Γfs_mask)
    Γd2_to_Γfs = findall(Γd2_to_Γfs_mask)
    Γfs_mid_to_Γfs = findall(!, Γd1_to_Γfs_mask .| Γd2_to_Γfs_mask)

    Γd1 = Triangulation(Γfs, Γd1_to_Γfs)
    Γd2 = Triangulation(Γfs, Γd2_to_Γfs)
    Γfs_mid = Triangulation(Γfs, Γfs_mid_to_Γfs)

    return (
        Ω=Ω,
        Γfs=Γfs,
        Γfs_mid=Γfs_mid,
        Γd1=Γd1,
        Γd2=Γd2,
        Γstr=Γstr,
        Λstr=Λstr,
        nΛstr=nΛstr,
        Γbot=Γbot,
        Γinlet=Γinlet,
        Γoutlet=Γoutlet,
    )
end

"Building measures for each region the 2D problem domain"      # Measure defines how to integrate over the domain and its boundaries
function build_measures_khabakhpasheva(reg; order::Int)
    degree  = 2*order  
    dΩ      = Measure(reg.Ω, degree)
    dΓfs    = Measure(reg.Γfs, degree)
    dΓfs_mid = Measure(reg.Γfs_mid, degree)
    dΓd1     = Measure(reg.Γd1, degree)
    dΓd2     = Measure(reg.Γd2, degree)
    dΓstr   = Measure(reg.Γstr, degree)
    dΓbot   = Measure(reg.Γbot, degree)
    dΓinlet = Measure(reg.Γinlet, degree)
    dΓoutlet= Measure(reg.Γoutlet, degree)
    dΛstr   = Measure(reg.Λstr, degree)

    return (
        dΩ=dΩ,
        dΓfs=dΓfs,
        dΓfs_mid=dΓfs_mid,
        dΓd1=dΓd1,
        dΓd2=dΓd2,
        dΓstr=dΓstr,
        dΓbot=dΓbot,
        dΓinlet=dΓinlet,
        dΓoutlet=dΓoutlet,
        dΛstr=dΛstr
    ) 
end