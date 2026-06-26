# domain.jl
#
# Builds the Gridap triangulations (regions) and integration measures
# for each case's problem domain, from a loaded mesh model. Case 1 and
# Case 3 split the free surface into inlet/outlet damping zones plus a
# central measurement zone, used to absorb outgoing waves in the
# frequency-domain formulation. Case 2's time-domain formulation has no
# damping zones; outgoing-wave reflection is instead limited by using a
# longer physical domain.

using GridapGmsh

"Load a Gmsh .msh file into a Gridap discrete model"
function build_model(filename)
    model = GmshDiscreteModel(filename)
    return model
end

"Building tagged regions for the case 1 problem domain"         # Defines where cells and faces are located in the domain, and how they are connected
function build_case1_regions(model, params)
    Ω = Triangulation(model)                                
    Γbot    = BoundaryTriangulation(model, tags=TAG_bot)
    Γinlet  = BoundaryTriangulation(model, tags=TAG_inlet)
    Γoutlet = BoundaryTriangulation(model, tags=TAG_outlet)
    Γfs  = BoundaryTriangulation(model, tags=TAG_fs)
    Γstr = BoundaryTriangulation(model, tags=TAG_str)

    Γ = BoundaryTriangulation(model)
    n = get_normal_vector(Γ)
    Λstr  = Skeleton(Γstr)
    nΛstr = get_normal_vector(Λstr)

    xΓfs = get_cell_coordinates(Γfs)

    function is_damping1(xs)
        npts = length(xs)
        xmid = sum(x[1] for x in xs) / npts
        return params.x₀ <= xmid <= params.xdᵢₙ
    end

    function is_damping2(xs)
        npts = length(xs)
        xmid = sum(x[1] for x in xs) / npts
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
        n=n
    )
end    # Returns a NamedTuple

"Building measures for each region of the case 1 problem domain"   # Measure defines how to integrate over the domain and its boundaries
function build_case1_measures(reg; order::Int)
    degree  = 2*order                                       # The accuracy of the integration
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
        dΛstr=dΛstr,
        dΓbot=dΓbot,
        dΓinlet=dΓinlet,
        dΓoutlet=dΓoutlet
    ) # Returns a NamedTuple (ordered list)
end

"Building tagged regions for the time-dependent case 2 problem domain"  # Defines where cells and faces are located in the domain, and how they are connected
function build_case2_transient_regions(model)
    Ω    = Triangulation(model)
    Ωp   = Triangulation(model, tags=TAG_porous)
    Γfs  = BoundaryTriangulation(model, tags=TAG_fs)
    Γstr = BoundaryTriangulation(model, tags=TAG_str)
    Λstr = SkeletonTriangulation(Γstr)

    nfs   = get_normal_vector(Γfs)
    nstr  = get_normal_vector(Γstr)
    nΛstr = get_normal_vector(Λstr)

    return (; Ω, Ωp, Γfs, Γstr, Λstr, nfs, nstr, nΛstr, model)
end

"Building measures for each region of the time-dependent case 2 problem domain"   # Measure defines how to integrate over the domain and its boundaries
function build_case2_transient_measures(reg; order::Int)
    degree = 2*order
    dΩ    = Measure(reg.Ω,    degree)
    dΩp   = Measure(reg.Ωp,   degree)
    dΓfs  = Measure(reg.Γfs,  degree)
    dΓstr = Measure(reg.Γstr, degree)
    dΛstr = Measure(reg.Λstr, degree)
    return (; dΩ, dΩp, dΓfs, dΓstr, dΛstr)
end

"Building tagged regions for the case 3 problem domain"           # Defines where cells and faces are located in the domain, and how they are connected
function build_case3_regions(model, params)
    Ω = Triangulation(model)
    Γbot    = BoundaryTriangulation(model, tags=TAG_bot)
    Γinlet  = BoundaryTriangulation(model, tags=TAG_inlet)
    Γoutlet = BoundaryTriangulation(model, tags=TAG_outlet)
    Γfs     = BoundaryTriangulation(model, tags=TAG_fs)
    Γpipe   = BoundaryTriangulation(model, tags=TAG_pipe)

    Γ = BoundaryTriangulation(model)
    n = get_normal_vector(Γ)
    nfs   = get_normal_vector(Γfs)
    npipe = get_normal_vector(Γpipe)

    xΓfs = get_cell_coordinates(Γfs)

    function is_damping1(xs)
        npts = length(xs)
        xmid = sum(x[1] for x in xs) / npts
        return params.x₀ <= xmid <= params.xdᵢₙ
    end

    function is_damping2(xs)
        npts = length(xs)
        xmid = sum(x[1] for x in xs) / npts
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
        Γpipe=Γpipe,
        Γbot=Γbot,
        Γinlet=Γinlet,
        Γoutlet=Γoutlet,
        n=n,
        nfs=nfs,
        npipe=npipe
    )
end

"Building measures for each region of the case 3 problem domain"   # Measure defines how to integrate over the domain and its boundaries
function build_case3_measures(reg; order::Int)
    degree   = 2*order                                       # The accuracy of the integration
    dΩ       = Measure(reg.Ω, degree)
    dΓfs     = Measure(reg.Γfs, degree)
    dΓfs_mid = Measure(reg.Γfs_mid, degree)
    dΓd1     = Measure(reg.Γd1, degree)
    dΓd2     = Measure(reg.Γd2, degree)
    dΓbot    = Measure(reg.Γbot, degree)
    dΓinlet  = Measure(reg.Γinlet, degree)
    dΓoutlet = Measure(reg.Γoutlet, degree)
    dΓpipe   = Measure(reg.Γpipe, degree)

    return (
        dΩ=dΩ,
        dΓfs=dΓfs,
        dΓfs_mid=dΓfs_mid,
        dΓd1=dΓd1,
        dΓd2=dΓd2,
        dΓbot=dΓbot,
        dΓinlet=dΓinlet,
        dΓoutlet=dΓoutlet,
        dΓpipe=dΓpipe
    ) # Returns a NamedTuple (ordered list)
end
