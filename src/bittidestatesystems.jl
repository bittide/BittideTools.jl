
module BittideStateSystems

using JuliaTools

# This is a dependency of Callisto, it does not depend on Callisto.

export OneEdgeOutputSystem, OutputSystem, Measurement, make_output_system

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

function StateSystems.next(K::OneEdgeOutputSystem, measurement)
    measurement.occupancies[K.portnum] - K.offset
end

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

########################################################################
# controller measurement

# IOM = instant of measurement
# IOFC = instant of frequency change

##############################################################################

##############################################################################
# We use this object to preallocate storage for the measured occupancies

mutable struct Measurement
    occupancies::Vector{Float64} # indexed by port number
    theta_at_ioc::Float64
    physical_time_at_iom::Float64
    previous_incoming_link_status::Vector{Int64}
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

function StateSystems.next(K::OutputSystem, measurement)
    sum(measurement.previous_incoming_link_status .*
        (measurement.occupancies - K.local_offsets))
end

end
