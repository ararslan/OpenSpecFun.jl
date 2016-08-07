__precompile__()

module OpenSpecFun

import Base.Math: @horner, @pg_horner

deps = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
if isfile(deps)
    include(deps)
else
    error("OpenSpecFun.jl not properly installed. Please run Pkg.build(\"OpenSpecFun\").")
end

for f in ["bessel", "erf", "gamma"]
    include("$f.jl")
end

end # module
