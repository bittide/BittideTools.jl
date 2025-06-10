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

export CalOpts, check_freq_is_positive, Checker, Error, errorat,
    get_beta_constant_ugn, get_constant_ugn, get_frequencies,
    get_latency, get_offset, initial_error,
    local_to_realtime, make_frequencies,
    realtime_to_local, sim_is_done

"""
    CalOpts structure
Holds configuration options for a Callisto simulation.

# Fields
- `graph`:                      The network topology. Defaults to a 3-node bidirectional triangle.
- `latency`:                    Link latencies. Can be a scalar (applied to all links) or a vector. Defaults to 20.
- `offset`:                     Initial buffer offsets/target occupancies. Can be a scalar or a vector. Defaults to 200.
- `tmax`:                       Maximum simulation time. Defaults to 1.5e5.
- `kp`:                         Proportional gain for the PI controller. Defaults to 8e-6.
- `ki`:                         Integral gain for the PI controller. Defaults to 0.0.
- `poll_period`:                Time between controller updates for each node. Defaults to 2000.
- `theta0`:                     Initial phase for each node. Can be a scalar or a vector. Defaults to 0.0.
- `errors`:                     Clock frequency errors for each node.
                                Can be `:random` (generates random errors), a vector of numerical error values,
                                or a vector of `SimCore.Error` objects. Defaults to `:random`.
- `base_freq`:                  Base frequency used by the controller for scaling the integral term. Defaults to 1.0.
- `stopper`:                    A custom stopping condition for the simulation. Defaults to `nothing` (no custom stopper).
- `check_freq_is_positive`:     If `true`, checks if node frequencies become negative during simulation. Defaults to `true`.
- `enable_logging`:             If `true`, enables logging of simulation data (e.g., buffer occupancies). Defaults to `true`.
- `controller`:                 Specifies the controller type for each node.
                                Can be `:pi` (default PI controller) or a vector of
                                custom controller objects. Defaults to `:pi`.
- `v`:                          Options specific to Callisto1 or Callisto2 features (e.g., link up/down behavior, prefill).


"""
Base.@kwdef mutable struct CalOpts
    graph = Topology.Graph(["triangle"]; bidirectional = true)
    latency = 5000
    offset = 50
    tmax = 1e9
    kp = 2e-8
    ki = 1e-15
    poll_period = 100000
    theta0 = 0.1
    errors = :random
    base_freq = 1.0
    stopper = nothing
    check_freq_is_positive = true
    enable_logging = true
    controller = :pi
    v = nothing
end
# CalOpts should be such that at any time before running callisto,
# you can change any of the entries of CalOpts. That is, the constructor
# should not enforce constraints between the fields. That is why
# we do not instantiate the errors with a random vector of length
# n, since if subsequently the graph is changed then the struct
# would become inconsistent.



# convenience
CalOpts(graph::Topology.Graph; kw...) = CalOpts(;graph, kw...)


function make_frequencies(seed, num_nodes)
    rng = Random.Xoshiro(seed)
    f = 1 .+  rand(rng, 1:10^5, num_nodes) *  1e-9
    return f
end

function get_frequencies(c::CalOpts)
    if c.errors == :random
        return make_frequencies(1, c.graph.n)
    end
    return c.errors
end

function get_latency(c::CalOpts)
    latency = [ati(c.latency, e) for e=1:c.graph.m]
    return latency
end

function get_offset(c::CalOpts)
    offset = [ati(c.offset, e) for e=1:c.graph.m]
    return offset
end

beta_internal(ugn, src, dst, latency, gear, t, theta) = ugn + floor(gear*theta[src](t - latency)) - floor(gear*theta[dst](t))

function get_beta_constant_ugn(c, e, t, theta)
    src = c.graph.edges[e].src
    dst = c.graph.edges[e].dst
    beta0 = get_offset(c)
    ugn = get_constant_ugn(c, beta0, errors)
    latency = get_latency(c)
    return beta_internal(ugn[e], src, dst, latency[e], 1, t, theta)
end


function mkugn_internal(beta0, gear, latency, theta0_at_src,  wm2_at_src, theta0_at_dst)
    return beta0 - floor(gear*(theta0_at_src - latency * wm2_at_src)) + floor(gear*theta0_at_dst)
end

# only works for Callisto v1 style UGNs
# which are constant over time.
#
# Since this may be called by other code (e.g. CallistoLinear)
# which does not have gears, we make the gears
# an option.
function get_constant_ugn(c::CalOpts, beta0, errors; gears = 1)
    m = c.graph.m
    n = c.graph.n
    theta0 = [ati(c.theta0, i) for i = 1:n]
    gears = [ati(gears, e) for e=1:m]
    beta0 = [ati(beta0, e) for e=1:m]
    latencies = [ati(c.latency, e) for e = 1:m]
    wm2 = [initial_error(errors[i]) for i=1:n]
    ugn = [mkugn_internal(
        beta0[e],
        gears[e],
        latencies[e],
        theta0[c.graph.edges[e].src],
        wm2[c.graph.edges[e].src],
        theta0[c.graph.edges[e].dst]) for e=1:m]
    return ugn
end








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