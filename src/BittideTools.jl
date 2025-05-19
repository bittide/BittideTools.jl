
module BittideTools

#
# This is a dependency of Callisto.
#

include("bittidestatesystems.jl")
using .BittideStateSystems

include("bittidematrices.jl")
using .BittideMatrices

function reexport(m)
    for a in names(m)
        eval(Expr(:export, a))
    end
end

reexport(BittideStateSystems)
reexport(BittideMatrices)



end
