# Wave Attenuation by Floating Ecosystems for Floating Photovoltaic Systems

Code accompanying the MSc thesis *Wave Attenuation by Floating Ecosystems for
Floating Photovoltaic Systems* (TU Delft, 2026), developed with industry
partner Solinoor.

The model evaluates three representations of a floating wave barrier at the
Willemspolder site:

- **Case 1** — thin elastic Euler-Bernoulli beam, potential flow, frequency
  domain, unmoored.
- **Case 2** — elastic beam with a submerged porous root zone, time domain,
  unmoored.
- **Case 3** — rigid hollow pipe, potential flow, frequency domain, moored
  and unmoored.

All three are implemented with [Gridap.jl](https://github.com/gridap/Gridap.jl).

Case 1's formulation follows the monolithic FEM approach of Colomés,
Verdugo and Akkerman (2023), coupling potential flow to an Euler-Bernoulli
beam via a continuous/discontinuous Galerkin discretisation. The reference
implementation is [`MonolithicFEMVLFS.jl`](https://github.com/oriolcg/MonolithicFEMVLFS.jl).
The Khabakhpasheva benchmark validation scripts in this repository
(`case1_khab_validation.jl`, `case2_khab_validation.jl`) reproduce the same
test case used there.

> Colomés, O., Verdugo, F., & Akkerman, I. (2023). A monolithic finite
> element formulation for the hydroelastic analysis of very large floating
> structures. *International Journal for Numerical Methods in Engineering*,
> 124(3), 714–751.

## Requirements

- Julia 1.12 (developed and tested with 1.12.5)
- A working Gmsh installation, via `GridapGmsh.jl`

From the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Each script also activates and instantiates the project itself, so this step
is mainly useful for confirming the environment resolves before running
anything.

## Structure

```
src/
  Thesis_FPVBarrier.jl        Main module
  core/                       Mesh loading, regions, measures, FE spaces,
                              boundary tags — shared across all three cases
  forms/                      Weak forms / residuals for each case, plus the
                              Khabakhpasheva benchmark validation operators

scripts/
  case1.jl                     Case 1 base-case solve (single frequency)
  case3.jl                     Case 3 base-case solve (single frequency)
  case2_sweep_L4.jl            Case 2 time-domain sweep, L_b = 4.6 m
  case2_sweep_L18.jl           Case 2 time-domain sweep, L_b = 18.4 m
  case1_khab_validation.jl     Case 1/3 validation against Khabakhpasheva et al.
  case2_khab_validation.jl     Case 2 validation against Khabakhpasheva et al.

postprocess/
  case1_postprocess.jl          Case 1 parameter sweeps and figures
  case3_postprocess.jl          Case 3 parameter sweeps and figures
  case2_time_postprocess.jl     Case 2 figures, reanalysed from saved JLD2 output
  nitsche_validation_td.jl      Nitsche kinematic-constraint check, Case 2
  spectral_weighting.jl         Spectrally-weighted load reduction (ΔE_L,spec)
  structural_impl.jl            Load and fatigue reduction factors, cycle counting
  wave_extraction.jl            Shared wave amplitude extraction helpers

meshes/
  case1_length_sweep/           L_b ∈ {4.6, 9.2, 18.4} m
  case2_length_sweep/           L_b ∈ {4.6, 18.4} m
  case3_diametersweep/          D ∈ {0.25, ..., 0.50} m
  validation/                   Khabakhpasheva benchmark meshes

results/
  data/ref_data/                 Reference data (Khabakhpasheva, Riyansyah)
                                  used by the validation scripts
```

`results/jld2/`, `results/plots/`, and `results/vtk/` are created by the
scripts themselves and are not tracked in the repository.

## Running

Each case follows the same pattern: a base script builds the model and runs
a single solve, and a postprocess script (run afterwards, in the same Julia
session) performs the parameter sweeps and produces the figures and tables
used in the thesis.

**Case 1:**
```julia
include("scripts/case1.jl")
include("postprocess/case1_postprocess.jl")
```

**Case 2** (each sweep script is fully self-contained, no postprocess
session dependency):
```
julia scripts/case2_sweep_L4.jl
julia scripts/case2_sweep_L18.jl
julia postprocess/case2_time_postprocess.jl
```

**Case 3:**
```julia
include("scripts/case3.jl")
include("postprocess/case3_postprocess.jl")
```

**Validation:**
```
julia scripts/case1_khab_validation.jl
julia scripts/case2_khab_validation.jl
julia postprocess/nitsche_validation_td.jl
```

**Ch. 6 implications for FPV systems** (run after Case 1 and Case 3's
postprocess scripts, since these load the saved sweep results):
```
julia postprocess/spectral_weighting.jl
julia postprocess/structural_impl.jl
```

## Notes

- Case 2's quantitative KT/KR extraction was attempted but is not used in
  the thesis: outlet-boundary reflections corrupt the downstream probe
  signal, making the decomposition unreliable. Case 2 results are reported
  qualitatively, via free-surface and beam-deflection profiles.
- The `*_khab_validation` scripts and `nitsche_validation_td.jl` validate
  the implementation against a published benchmark and a numerical
  consistency check respectively, independent of the Willemspolder site
  results.
