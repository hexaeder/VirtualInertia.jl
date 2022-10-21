var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = VirtualInertia","category":"page"},{"location":"#VirtualInertia","page":"Home","title":"VirtualInertia","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for VirtualInertia.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [VirtualInertia]","category":"page"},{"location":"#NetworkDynamics.ODEEdge","page":"Home","title":"NetworkDynamics.ODEEdge","text":"Conventions: (Anti-)symmetric lines:\n\ninputs:   u_src_r, u_src_i, u_dst_r, u_dst_i\noutputs:  i_r, i_i\ncurrent direction is defined from src to dst\ndst node will receive -i_r, -i_i\n\nAsymmetric lines:\n\ninputs:   u_src_r, u_src_i, u_dst_r, u_dst_i\noutputs:  i_src_r, i_src_i, i_dst_r, i_dst_i\nnot yet implemented. might by tricky with fidutial...\n\n\n\n\n\n","category":"type"},{"location":"#MetaGraphs.get_prop-Tuple{Any, Union{Graphs.AbstractEdgeIter, AbstractUnitRange, Vector}, Symbol}","page":"Home","title":"MetaGraphs.get_prop","text":"get_prop(g, keys::Iterable, prop::Symbol)\n\nGet same property prop with different values vals for differet keys.\n\n\n\n\n\n","category":"method"},{"location":"#MetaGraphs.set_prop!-Tuple{Any, Union{Graphs.AbstractEdgeIter, AbstractUnitRange, Vector}, Symbol, Vector}","page":"Home","title":"MetaGraphs.set_prop!","text":"set_prop!(g, keys::Iterable, prop::Symbol, vals::Iterable)\n\nSet same property prop with different values vals for differet identifiers keys.\n\n\n\n\n\n","category":"method"},{"location":"#VirtualInertia.DroopControl-Tuple{}","page":"Home","title":"VirtualInertia.DroopControl","text":"DroopControl(; params...)\n\nReturn block for droop control outer ocntrol.\n\n\n\n\n\n","category":"method"},{"location":"#VirtualInertia.PT1Source-Tuple{}","page":"Home","title":"VirtualInertia.PT1Source","text":"PT1Source(; params...)\n\nCreate Schiffer Voltage source which follows angle directly but amplitude with lag.\n\n\n\n\n\n","category":"method"},{"location":"#VirtualInertia.PerfectSource-Tuple{}","page":"Home","title":"VirtualInertia.PerfectSource","text":"PerfectSource()\n\nPerfect Voltage source which follows the reference directly.\n\n\n\n\n\n","category":"method"},{"location":"#VirtualInertia._group_states-Tuple{Any}","page":"Home","title":"VirtualInertia._group_states","text":"Find complex states with _r and _i suffix\n\n\n\n\n\n","category":"method"},{"location":"#VirtualInertia.bustype-Tuple{Int64}","page":"Home","title":"VirtualInertia.bustype","text":"bustype(i::Int)\nbustype(s::String)\nbustype(s::Symbol)\n\nHelper to convert the different representations of bus types into eachother. Base representation is Symbol:\n\n:load   == bustype(1) == bustype(\"PQ\") == bustype(:PQ)\n:gen    == bustype(2) == bustype(\"PV\") == bustype(:PV)\n:syncon == bustype(3) == bustype(\"Ref\") == bustype(:Ref)\n\n\n\n\n\n","category":"method"},{"location":"#VirtualInertia.describe_edges-Tuple{MetaGraphs.MetaGraph}","page":"Home","title":"VirtualInertia.describe_edges","text":"describe_nodes(g::MetaGraph; firstcols=Vector{String}())\n\nReturns DataFrame with all the edge meta data.\n\n\n\n\n\n","category":"method"},{"location":"#VirtualInertia.describe_nodes-Tuple{MetaGraphs.MetaGraph}","page":"Home","title":"VirtualInertia.describe_nodes","text":"describe_nodes(g::MetaGraph; firstcols=Vector{String}())\n\nReturns DataFrame with all the node meta data.\n\n\n\n\n\n","category":"method"},{"location":"#VirtualInertia.import_system-Tuple{Symbol}","page":"Home","title":"VirtualInertia.import_system","text":"import_system(sym::Symbol; kwargs...)::MetaGraph\n\nMain entry point to load the systems. New systems should overload this function. Known implementations\n\nimport_system(:rtsgmlc): loads the GMLC update for the rts96\n\nThose functions return a MetaGraph which has properties attached to the Nodes/Edges.\n\nGraph properties:\n\nPbase : Base power for PU\n\noptional:\n\nNodeProps : tuple of node property names which should appear first in describe functions\nEdgeProps : tuple of edge property names which should appear first in describe functions\n\nNode properties:\n\nP : active power in PU\nQ : reactive power in PU\nVbase : base voltage in kV\nVm : voltage magnitude in PU\n\noptional:\n\ninertia : inertia in MJ/MW (see H in Wikipedia page)\ndamping : damping factor γ for swing equations\ntimeconstant : time cosntant τ for dynamic loads\nx, y : position of Bus for plotting purposes\nVa : voltage angle in rad\n\nEdge properties:\n\nrating : short term emergency rating in PU\nR : resistance in PU\nX : reactance in PU\n\n\n\n\n\n","category":"method"},{"location":"#VirtualInertia.import_system-Tuple{Val{:rtsgmlc}}","page":"Home","title":"VirtualInertia.import_system","text":"import_system(:rtsgmlc)\n\nImport the RTS-GMLC system as a MetaGraph.\n\n\n\n\n\n","category":"method"},{"location":"#VirtualInertia.subscript-Union{Tuple{T}, Tuple{T, Int64}} where T<:Union{AbstractString, Symbol}","page":"Home","title":"VirtualInertia.subscript","text":"subscript(s, i)\n\nAppend symbol or string s with a integer subscript.\n\n\n\n\n\n","category":"method"}]
}
