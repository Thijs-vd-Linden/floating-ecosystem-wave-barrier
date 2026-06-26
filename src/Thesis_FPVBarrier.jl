module Thesis_FPVBarrier

using LinearAlgebra
using SparseArrays
using Gridap
using Plots

# Include source files
include("core/parameters.jl")
include("core/tags.jl")
include("core/domain.jl")
include("core/spaces.jl")


include("forms/case1_form.jl")
include("forms/case2_transient_form.jl")
include("forms/case3_form.jl")

include("forms/case1_khab_validation_operator.jl") 
include("forms/case2_validation_operator.jl")
include("core/khabakhpasheva_domain.jl")    

include("../postprocess/wave_extraction.jl")
end 