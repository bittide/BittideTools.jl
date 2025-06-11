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

export check_freq_is_positive, Checker, Error, errorat,
    get_beta_constant_ugn, get_constant_ugn, get_frequencies,
    get_latency, get_offset, initial_error,
    local_to_realtime, make_frequencies,
    realtime_to_local, sim_is_done





#CalOpts(graph::Topology.Graph; kw...) =  CalParams(; graph, kw...)



# Base.@kwdef mutable struct ComputedParams
#     errors
#     latency
#     offset
#     theta0
#     c1_beta0
#     c1_gears
#     c1_control_delay
# end


#
# Have two separate CalOpts structs. Not just one.
# And a third one for Odes. Etc.
# and converters between them.
#
#



#
#
# function make_link(e)
#     return Link(e, g.edges[e].src, g.edges[e].dst,
#                 make_ugn(beta0[e], gears[e], latency[e], theta0[g.edges[e].src],
#                     wm2[g.edges[e].src], theta0[g.edges[e].dst]),
#                 latency[e],
#                 gears[e],
#                 beta0[e],
#                 measurement_offsets[e])
# end
#
# links = [make_link(e) for e = 1:g.m]

# # convenience function
# function Base.getproperty(cp::CalParams, sym::Symbol)
#     if sym == :links
#         return cp.model.links
#     elseif sym == :errors
#         if !isnothing(cp.model)
#             return cp.model.errors
#         end
#         return getfield(cp, :errors)
#     end
#     return getfield(cp, sym)
# end

# function make_frequencies(seed, num_nodes)
#     rng = Random.Xoshiro(seed)
#     f = 1 .+  rand(rng, 1:10^5, num_nodes) *  1e-9
#     return f
# end

# function get_frequencies(c::CalOpts)
#     if c.errors == :random
#         return make_frequencies(1, c.graph.n)
#     end
#     return c.errors
# end

# function get_latency(c::CalOpts)
#     latency = [ati(c.latency, e) for e=1:c.graph.m]
#     return latency
# end

# function get_offset(c::CalOpts)
#     offset = [ati(c.offset, e) for e=1:c.graph.m]
#     return offset
# end

# beta_internal(ugn, src, dst, latency, gear, t, theta) = ugn + floor(gear*theta[src](t - latency)) - floor(gear*theta[dst](t))

# function get_beta_constant_ugn(c, e, t, theta)
#     src = c.graph.edges[e].src
#     dst = c.graph.edges[e].dst
#     beta0 = get_offset(c)
#     ugn = get_constant_ugn(c, beta0, errors)
#     latency = get_latency(c)
#     return beta_internal(ugn[e], src, dst, latency[e], 1, t, theta)
# end


# function mkugn_internal(beta0, gear, latency, theta0_at_src,  wm2_at_src, theta0_at_dst)
#     return beta0 - floor(gear*(theta0_at_src - latency * wm2_at_src)) + floor(gear*theta0_at_dst)
# end

# only works for Callisto v1 style UGNs
# which are constant over time.
#
# Since this may be called by other code (e.g. CallistoLinear)
# which does not have gears, we make the gears
# # an option.
# function get_constant_ugn(c::CalOpts, beta0, errors; gears = 1)
#     m = c.graph.m
#     n = c.graph.n
#     theta0 = [ati(c.theta0, i) for i = 1:n]
#     gears = [ati(gears, e) for e=1:m]
#     beta0 = [ati(beta0, e) for e=1:m]
#     latencies = [ati(c.latency, e) for e = 1:m]
#     wm2 = [initial_error(errors[i]) for i=1:n]
#     ugn = [mkugn_internal(
#         beta0[e],
#         gears[e],
#         latencies[e],
#         theta0[c.graph.edges[e].src],
#         wm2[c.graph.edges[e].src],
#         theta0[c.graph.edges[e].dst]) for e=1:m]
#     return ugn
# end








################################################################################

struct Checker end

@inline check_freq_is_positive(ch::Nothing, c, errors, i, s) = return

function check_freq_is_positive(ch::Checker, c, errors, i, s)
    if c + errorat(errors[i], s) <= 0
        println("Frequency at node $i has become negative: omega = ", c + errors[i](s))
        println("time s = ", s)
    end
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