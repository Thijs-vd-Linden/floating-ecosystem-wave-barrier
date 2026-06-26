# case3.jl
#
# Single base-case solve for Case 3: builds the model for the base
# diameter (D = 0.45 m), assembles the monolithic block system coupling
# the fluid finite-element field to the rigid-body heave/pitch degrees
# of freedom, and solves it directly. The block solution vector is then
# split into its FE part (potential, free-surface elevation) and rigid-
# body part (heave, pitch), and the pipe boundary displacement is
# reconstructed from the two rigid-body modes for VTK output.
#
# Run case3_postprocess.jl afterwards, in the same session, for the
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

vtk_dir = "results/vtk/case3"
mkpath(vtk_dir)
plot_dir = "results/plots/case3"
mkpath(plot_dir)
filename3 = joinpath(vtk_dir, "case3")

# settings
order = 2

# build mesh/model
model3   = Thesis_FPVBarrier.build_model("meshes/case3_diametersweep/case3_D_45.msh")
params3 = Thesis_FPVBarrier.build_case3_params(model3)                          # Parameters
reg3     = Thesis_FPVBarrier.build_case3_regions(model3, params3)
meas3    = Thesis_FPVBarrier.build_case3_measures(reg3; order=order)
sp3      = Thesis_FPVBarrier.build_case3_spaces(reg3, params3; order=order)
sys3     = Thesis_FPVBarrier.build_case3_operator(sp3, reg3, meas3, order, params3)

x = sys3.A \ sys3.b           # solves Ax=b
nfe = num_free_dofs(sp3.X)   # Extract the fluid and structure parts of the solution vector x
xfe = x[1:nfe]              # FEM solution for the fluid (ϕ, κ)
ξ   = x[nfe+1:end]          # Rigid body DOFs solution

nϕ = num_free_dofs(sp3.U_Ω)
nκ = num_free_dofs(sp3.U_κ)

println("nϕ = ", nϕ)
println("nκ = ", nκ)
println("length(xfe) = ", length(xfe))

xϕ = xfe[1:nϕ];
xκ = xfe[nϕ+1:nϕ+nκ];

ϕh3 = FEFunction(sp3.U_Ω, xϕ);
κh3 = FEFunction(sp3.U_κ, xκ);

# Postprocess
xc = params3.xc
zc = params3.zc
npipe = reg3.npipe;

ψ₁(x) = VectorValue(0.0, 1.0);
ψ₂(x) = VectorValue(-(x[2] - zc), x[1] - xc);
ψ₂_field = CellField(ψ₂, reg3.Γpipe);

H1 = npipe ⋅ ψ₁;
H2 = npipe ⋅ ψ₂_field;

ηph3 = ξ[1] * H1 + ξ[2] * H2;

# write results to VTK files for visualization in Paraview
writevtk(reg3.Ω, filename3 * "_fluid", cellfields = ["phi_re" => real(ϕh3), "phi_im" => imag(ϕh3), "phi_abs" => abs(ϕh3) ]);
writevtk(reg3.Γfs, filename3 * "_free_surface", cellfields = ["kappa_re"=>real(κh3), "kappa_im"=>imag(κh3), "kappa_abs"=>abs(κh3)]);
writevtk(reg3.Γpipe, filename3 * "_pipe", cellfields = ["eta_p_re"  => real(ηph3), "eta_p_im"  => imag(ηph3), "eta_p_abs" => abs(ηph3)]);

# DOFs
println("ξ_heave = ", ξ[1])
println("ξ_pitch = ", ξ[2])

