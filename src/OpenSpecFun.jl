__precompile__()

module OpenSpecFun

deps = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
if isfile(deps)
    include(deps)
else
    error("OpenSpecFun.jl not properly installed. Please run Pkg.build(\"OpenSpecFun\").")
end

import Base.Math: @horner, @pg_horner

# export
#     # bessel
#     airy,
#     airyai,
#     airyaiprime,
#     airybi,
#     airybiprime,
#     airyprime,
#     airyx,
#     besselh,
#     besselhx,
#     besseli,
#     besselix,
#     besselj,
#     besselj0,
#     besselj1,
#     besseljx,
#     besselk,
#     besselkx,
#     bessely,
#     bessely0,
#     bessely1,
#     besselyx,
#     hankelh1,
#     hankelh1x,
#     hankelh2,
#     hankelh2x,
#     # erf
#     erf,
#     erfc,
#     erfcinv,
#     erfcx,
#     erfi,
#     erfinv

const libm = Base.libm_name
const openspecfun = "libopenspecfun"

for f in ["bessel", "erf"]
    include("$f.jl")
end

end # module
