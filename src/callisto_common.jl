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

export CalOpts, make_frequencies


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
    latency =20
    offset = 200
    tmax = 1.5e5
    kp = 8e-6
    ki = 0.0
    poll_period = 2000
    theta0 = 0.0
    errors = :random
    base_freq = 1.0
    stopper = nothing
    check_freq_is_positive = true
    enable_logging = true
    controller = :pi
    v = nothing
end



function make_frequencies(seed, num_nodes)
    rng = Random.Xoshiro(seed)
    f = 1 .+  rand(rng, 1:10^5, num_nodes) *  1e-9
    return f
end






end