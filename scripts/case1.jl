# case1.jl
#
# Single base-case solve for Case 1: builds the model from one of the
# three length-sweep meshes (L_4 active by default), solves the
# frequency-domain operator at the design frequency, and writes the
# velocity potential, free-surface elevation, and beam displacement
# fields to VTK for inspection in ParaView.
#
# Run case1_postprocess.jl afterwards, in the same session, for the
# parameter sweeps and figures used in the thesis.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..")) 
Pkg.instantiate()

include("../src/Thesis_FPVBarrier.jl")
using .Thesis_FPVBarrier
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using CSV
using Plots

vtk_dir = "results/vtk/case1"
mkpath(vtk_dir)
plot_dir = "results/plots/case1"
mkpath(plot_dir)
filename1 = joinpath(vtk_dir, "case1")

# settings
order = 2

# build mesh/model, uncomment the mesh you want to use for the simulation
model1 = Thesis_FPVBarrier.build_model("meshes/case1_length_sweep/case1_L_4.msh")
# model1 = Thesis_FPVBarrier.build_model("meshes/case1_length_sweep/case1_L_9.msh")
# model1 = Thesis_FPVBarrier.build_model("meshes/case1_length_sweep/case1_L_18.msh")

# Parameters
params1 = Thesis_FPVBarrier.build_case1_params(model1)
reg1   = Thesis_FPVBarrier.build_case1_regions(model1, params1)
meas1  = Thesis_FPVBarrier.build_case1_measures(reg1; order=order)
sp1    = Thesis_FPVBarrier.build_case1_spaces(reg1; order=order)
@time op1 = Thesis_FPVBarrier.build_case1_operator(sp1, reg1, meas1, order, params1)

(ϕh1, κh1, ηh1) = solve(op1);

writevtk(reg1.Ω, filename1 * "_fluid", cellfields = ["phi_re" => real(ϕh1), "phi_im" => imag(ϕh1), "phi_abs" => abs(ϕh1) ])
writevtk(reg1.Γfs, filename1 * "_free_surface", cellfields = ["kappa_re"=>real(κh1), "kappa_im"=>imag(κh1), "kappa_abs"=>abs(κh1)])
writevtk(reg1.Γstr, filename1 * "_structure", cellfields = ["eta_re"=>real(ηh1), "eta_im"=>imag(ηh1), "eta_abs"=>abs(ηh1)])


