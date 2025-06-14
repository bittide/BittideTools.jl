# Copyright 2023 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

module CallistoCommon

import Topology
using Random
using JuliaTools
using Piecewise

export check_freq_is_positive, check_initial, Checker, Error, errorat, get_beta_constant_ugn,
       get_constant_ugn, get_frequencies, get_latency, initial_error,
       local_to_realtime, make_frequencies, realtime_to_local, sim_is_done


function check_initial(p, w, n, tmax)
    avgfreq = sum(w)/n
    if avgfreq*tmax/p < 10
        println("Warning; fewer than 10 steps per node in simulation")
    end
    relmadf = sum(abs.(w .- avgfreq))/avgfreq
    if relmadf < 1e-9
        println("Warning: frequencies are very close together. Expect strange behavior.")
    end
end


function make_frequencies(seed, num_nodes)
     rng = Random.Xoshiro(seed)
     f = 1 .+  rand(rng, 1:10^5, num_nodes) *  1e-9
     return f
end





################################################################################

struct Checker end

@inline check_freq_is_positive(ch::Nothing, c, errors, i, s) = return true

function check_freq_is_positive(ch::Checker, c, errors, i, s)
    if c + errorat(errors[i], s) <= 0
        println("Frequency at node $i has become negative: omega = ", c + errorat(errors[i],s))
        println("time s = ", s)
        return false
    end
    return true
end


sim_is_done(stopper::Nothing, errors, i, p, c, s) = false

################################################################################
# errors


abstract type Error end

struct PwlError <: Error
    pwl::PiecewiseLinear
end

initial_error(e::Error) = e(0)
initial_error(e::Number) = e
initial_error(e::Vector{T}) where {T <: Error} = e(0)
initial_error(e::Vector{Float64}) = e
errorat(e::Number, t) = e
errorat(e::PwlError, t) = e.pwl(t)
errorat(e::Vector{T}, t) where {T <: Error}  =  [a(t) for a in e]
errorat(e::Vector{T}, t) where {T <: Number}  = e
Error(x::Number) = x
Error(p::PiecewiseLinear) = PwlError(p)
(e::PwlError)(t) = e.pwl(t)
(e::Vector{T})(t) where {T <: Error} = [a(t) for a in e]




"""
    local_to_realtime(e, p, c, s, wmin) -> Float64

Convert a local time interval `p` (in local clock ticks) to its equivalent duration `ds`
in wall-clock time.

This function solves the integral equation for ds

    int_s^{s+ds} (c + e(t)) dt = p

where `e(t)` is the clock frequency error at time `t` and `c` is a
constant frequency correction, and `s` is the starting wall-clock time.

# Arguments
- `e`: The clock error model. Can be a `Number` (for constant error) or a `PwlError`
       (for piecewise linear error).
- `p::Real`: The local time interval (duration in local ticks) to convert.
- `c::Real`: The constant frequency correction applied to the clock.
- `s::Real`: The starting wall-clock time.
- `wmin::Real`: A minimum effective frequency, used as a lower bound to prevent
                division by zero or excessively large `ds` calculations

# Returns
- `Float64`: The duration `ds` in wall-clock time.
"""
function local_to_realtime(e::PwlError, p, c, s, wmin)
    intfreq(s1, s2, c) = definite_integral(c + e.pwl, s1, s2)
    dt = bisection(s2 -> intfreq(s, s2, c) - p,  s, s + p/wmin) - s
    return dt
end
local_to_realtime(e::Number, p, c, s, wmin) = p / (c + e)

# convert a duration ds in localticks to a duration in wallclock
realtime_to_local(e::Number, c, ds) = ds * (c + e)








end