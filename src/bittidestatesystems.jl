
module BittideStateSystems

# Does not depend on Callisto.


export OneEdgeOutputSystem

using Topology
using StateSystems

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
OneEdgeOutputSystem(edgeid, c) = OneEdgeOutputSystem(c.graph, c.links[edgeid].offset, edgeid)




end

