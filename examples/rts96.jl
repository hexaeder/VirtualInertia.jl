using VirtualInertia
using MetaGraphs
using BlockSystems
using Unitful

network = import_system(:rtsgmlc)
nodes = describe_nodes(network)

sum(skipmissing(nodes.P_inj))
sum(nodes.P_load)
nodes.P_load

describe_edges(network)

# line model
piline = BSPiLine(ω0=2π*50, p_order=[:R, :L, :C])

# line parameters
R = ustrip.(u"pu", get_prop(network, edges(network), :R))
L = ustrip.(u"pu", get_prop(network, edges(network), :X))
# TODO: is this correct? That C = B? Maybe 1/B?
C = ustrip.(u"pu", get_prop(network, edges(network), :B))
lineP = zip(R, L, C) |> collect


unique(describe_nodes(network).type)
nodes = describe_nodes(network)
nodes[nodes.type .== :syncon, :]
nodes[nodes.type .== :load, :]
