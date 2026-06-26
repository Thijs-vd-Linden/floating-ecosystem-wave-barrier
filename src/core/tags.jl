# tags.jl
#
# Named constants for the mesh boundary and subdomain tags used across
# all three cases. These correspond to the physical group names defined
# in the .geo/.msh mesh files.

"Free Surface boundary"
const TAG_fs = "FreeSurface"

"Bottom boundary"
const TAG_bot = "Bottom"

"Inlet boundary (left)"
const TAG_inlet = "Inlet"

"Outlet boundary (right)"
const TAG_outlet = "Outlet"

"Fluid-structure boundary"
const TAG_str    = "Structure"

"Main fluid domain"
const TAG_fluid = "Fluid"

"Porous subdomain"
const TAG_porous = "Porous"

"Pipe structure"
const TAG_pipe = "Pipe"

"Free surface left edge (damping)"
const TAG_fsLeft = "FreeSurfaceLeft"

"Free surface right edge (damping)"
const TAG_fsRight = "FreeSurfaceRight"