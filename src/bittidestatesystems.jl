
module BittideStateSystems

using JuliaTools
import ..CallistoCommon: get_offset

# This is a dependency of Callisto, it does not depend on Callisto.


export OneEdgeOutputSystem, OutputSystem, Measurement, make_output_system, composed_controller

using Topology
import StateSystems: StateSystems, PIStateSystem, StateSystem, compose, next

mutable struct OneEdgeOutputSystem <: StateSystems.StateSystem
    portnum::Int64
    offset::Float64
    # This inner constructor only accepts arguments
    # of the specified types, so we can have an outer constructor
    # which takes more generic types
    OneEdgeOutputSystem(p::Int64, o::Number) = new(p, o)
end

StateSystems.next(K::OneEdgeOutputSystem, measurement) = measurement.occupancies[K.portnum] - K.offset


"""
    OneEdgeOutputSystem(graph, offset, edgeid)

Return a StateSystem which takes the measurement from Callisto
and returns the occupancy-offset of the edge. May only
be used as the controller at the destination node of the edge.
"""
function OneEdgeOutputSystem(graph::Topology.Graph, offset, edgeid)
    nodeid = graph.edges[edgeid].dst
    portnum = inportid_from_edgeid(graph::Graph, nodeid, edgeid)
    return OneEdgeOutputSystem(portnum, offset)
end

"""
    OneEdgeOutputSystem(edgeid, c)

Return a StateSystem which takes the measurement from Callisto
and returns the occupancy-offset of the edge. May only
be used as the controller at the destination node of the edge.
"""
OneEdgeOutputSystem(edgeid, c) = OneEdgeOutputSystem(c.graph, c.model.links[edgeid].offset, edgeid)


########################################################################
# controller measurement


# IOM = instant of measurement
# IOFC = instant of frequency change

##############################################################################




##############################################################################
# We use this object to preallocate storage for the measured occupancies

mutable struct Measurement
    occupancies::Vector{Float64} # indexed by port number
    theta_at_iofc::Float64
    physical_time_at_iom::Float64
    incoming_link_status::Vector{Int64}
end

"""
    Measurement(x)

Constructs a measurement object with occupancies defined by x.
Sets the initial link statuses to up.
"""
Measurement(x) = Measurement(x, 0, 0, ones(Int64, length(x)))

# This is a static system that converts the list
# of buffer occupancies into the required scalar for a P or PI controller
mutable struct OutputSystem <: StateSystem
    local_offsets::Vector{Float64}
end


StateSystems.next(K::OutputSystem, measurement) = sum(measurement.incoming_link_status .* (measurement.occupancies - K.local_offsets))

function make_output_system(i, offset, graph)
    local_offsets = inport_indexed_from_edge_indexed(graph, offset, i)
    return OutputSystem(local_offsets)
end

function composed_controller(i, graph, kp, ki, poll_period, base_freq, output_system)
    K1 = PIStateSystem(kp, ki * poll_period / base_freq)
    K2 = output_system
    return compose(K1, K2)
end

# # don't scale offsets by gears. instead
# # let the user deal with that.
# function make_local_offsets(i, c)
#     #offset = [c.model.links[e].offset for e = 1:c.graph.m]
#     offset = [ati(c.offset, e) for e = 1:c.graph.m]
#     local_offset = inport_indexed_from_edge_indexed(c.graph, offset, i)
#     return local_offset
# end


# c.computed.output_systems[i] = OutputSystem(make_local_offsets(i, c))
# default_controller(i, c) = default_controller(i, c.graph, c.kp, c.ki, c.poll_period, c.base_freq, make_local_offsets(i, c))



end

