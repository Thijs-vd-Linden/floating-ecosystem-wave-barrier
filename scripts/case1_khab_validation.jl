# case1_khab_validation.jl
#
# Validates Case 1 against the Khabakhpasheva et al. benchmark: builds
# the validation mesh and a spatially-varying EI profile matching the
# published reference geometry, solves the frequency-domain operator,
# and plots the resulting beam deflection against both the
# Khabakhpasheva and Riyansyah reference data.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..")) 
Pkg.instantiate()

include("../src/Thesis_FPVBarrier.jl")
using .Thesis_FPVBarrier
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using LaTeXStrings
using CSV
using Plots

vtk_dir = "results/vtk/validation/case1"
mkpath(vtk_dir)
plot_dir = "results/plots/validation/case1"
mkpath(plot_dir)
filename = joinpath(vtk_dir, "case1_khabakhpasheva")

# settings
order = 2

# Parameters
params = Thesis_FPVBarrier.build_khabakhpasheva_params()

# build mesh/model
model = Thesis_FPVBarrier.build_khabakhpasheva_model("meshes/validation/khabakhpasheva_domain.msh")
reg   = Thesis_FPVBarrier.build_regions_khabakhpasheva(model, params)
meas  = Thesis_FPVBarrier.build_measures_khabakhpasheva(reg; order=order)
sp    = Thesis_FPVBarrier.build_case1_spaces(reg; order=order)
op = Thesis_FPVBarrier.build_khabakhpasheva_operator(sp, reg, meas, order, params)

(ϕh, κh, ηh) = solve(op);

writevtk(reg.Ω, filename * "_fluid", cellfields = ["phi_re" => real(ϕh), "phi_im" => imag(ϕh), "phi_abs" => abs(ϕh) ])
writevtk(reg.Γfs, filename * "_free_surface", cellfields = ["kappa_re"=>real(κh), "kappa_im"=>imag(κh,), "kappa_abs"=>abs(κh)])
writevtk(reg.Γstr, filename * "_structure", cellfields = ["eta_re"=>real(ηh), "eta_im"=>imag(ηh), "eta_abs"=>abs(ηh)])

# Postprocess
xy_cp = get_cell_points(get_fe_dof_basis(sp.V_η)).cell_phys_point       # physical coordinates of the beam degrees of freedom
x_cp = [[xy_ij[1] for xy_ij in xy_i] for xy_i in xy_cp]                 # extract x-coordinates of the beam dofs
η_cdv = get_cell_dof_values(ηh)                                         # extract the corresponding η values at the beam dofs

# Flatten the per-cell lists so that we sort the entire beam once rather than
# basing the permutation on the first cell only, which can scramble later
# entries when multiple cells are present.
flat_x = vcat(x_cp...)
flat_η = vcat(η_cdv...)

p = sortperm(flat_x)                                                   # permutation that sorts the full x vector

xs = [(x - params.xb₀) / params.Lb for x in flat_x[p]]                # normalize the sorted coordinates as x/L
η_rel_xs = [abs(η) / params.η₀ for η in flat_η[p]]                    # relative beam displacement |η|/η₀, keeping order

# Reference data
ref_data1 = CSV.File("results/data/ref_data/Khabakhpasheva_without_joint.csv"; header=false)
ref_data2 = CSV.File("results/data/ref_data/Riyansyah_without_joint.csv"; header=false)

# Plot
plt = plot(
    xs,
    η_rel_xs,
    xlims=(0,1),
    lw=2,
    label="Case 1",
    xlabel=L"x/L_b",
    ylabel=L"|\eta|/\eta_0",
    legend=:topright,
    framestyle=:box,
    grid=true
)

plot!(
    plt,
    ref_data1.Column1,
    ref_data1.Column2,
    seriestype=:scatter,
    marker=:circle,
    ms=4,
    label="Khabakhpasheva et al."
)

plot!(
    plt,
    ref_data2.Column1,
    ref_data2.Column2,
    linestyle=:dash,
    lw=2,
    label="Riyansyah et al."
)

display(plt)
savefig(plt, joinpath(plot_dir, "case1_khabakhpasheva.png"))