
module BittideMatrices

using Topology
using LinearAlgebra



"""
   is_rate_matrix(Q)

Return true if Q is a rate matrix. That is, the offdiagonal elements are nonnegative
and the rows sum to zero.
"""
function is_rate_matrix(Q)
    n = size(Q,1)
    rowsums = sum(Q, dims = 2)
    if minimum([Q[i,j] for i=1:n, j=1:n if i != j]) < 0
        return false
    end
    maxelt  = maximum(abs.(Q))
    maxrowsums = maximum(abs.(rowsums))
    if maxrowsums / maxelt  > 1e-10
        return false
    end
    return true
end


# The matrix Q = DB' is a rate matrix. 
# The bittide dynamics is approximately \dot\theta = k Q \theta + v
#
function eigendecompose_metzler(A)
    @assert is_metzler(A) "Matrix is not Metzler"
    # given A Metzler
    n = size(A,1)
    F = eigen(A)
    p = sortperm(F.values, by=abs)
    T = F.vectors[:,p]
    c = real(T[1,1])
    T = T/c
    L = diagm(F.values[p])

    # Z = e^At is a stochastic matrix. Z\one = \one
    # e^At tends to 1 z^T
    z = real.(inv(T)[1,:])

    #@assert maximum(abs.(T*L*inv(T) - A)) < 1e-10
    #@assert maximum(abs.(exp(A*1e9) - ones(n)*z')) < 1e-12
    return T, L, z
end


function specinv(L)
    T, DD, z = eigendecompose_unsym(L)
    Di = inv(DD); Di[1,1]=0
    Li = real(T * Di * inv(T))
    return Li
end


# need stopper

end

