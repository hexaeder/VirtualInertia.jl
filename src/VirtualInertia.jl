module VirtualInertia

using NetworkDynamics
using BlockSystems
using ModelingToolkit
using ModelingToolkit: getname

include("utils.jl")

struct CallableBlockWrapper{F,G}
    f_ip::F
    g_oop::G
    states::Vector{Symbol}
    inputs::Vector{Symbol}
    rem_states::Vector{Symbol}
    params::Vector{Symbol}
end
function CallableBlockWrapper(gen)
    f_ip = gen.f_ip
    g_oop = gen.g_oop
    states = map(s->s.name, gen.states)
    inputs = map(s->s.name, gen.inputs)
    rem_states = map(s->s.name, gen.rem_states)
    params = map(s->s.name, gen.params)
    CallableBlockWrapper{typeof(f_ip), typeof(g_oop)}(f_ip,g_oop,states,inputs,rem_states,params)
end
(cbw::CallableBlockWrapper)(du,u,edges,p,t) = cbw.f_ip(du,u,flowsum(edges),p,t)

function NetworkDynamics.ODEVertex(inner::IOBlock, outer::IOBlock, p_order=[])
    cl = Inverter(inner, outer)
    open_p = ModelingToolkit.getname.(cl.iparams)
    if isempty(p_order) && !isempty(open_p)
        if length(open_p) > 1
            @warn "Model has open parameters $open_p. ODEVertex will use that order. You may specify manuallly using the `p_order` kw."
        end
        p_order = open_p
    end
    ODEVertex(cl, p_order)
end

function NetworkDynamics.ODEVertex(iob::IOBlock, p_order=[])
    @assert length(p_order) == length(iob.iparams) "Provide order of all iparams: $(getname.(iob.iparams))"

    if fulfills(iob, BlockSpec([:i_r, :i_i], [:u_r, :u_i]))
        # normal node
        gen = generate_io_function(iob;
                                   f_states=[iob.u_r, iob.u_i],
                                   f_inputs=[iob.i_r, iob.i_i],
                                   f_params=p_order,
                                   type=:ode,
                                   warn=false)
        f = CallableBlockWrapper(gen)
    elseif fulfills(iob, BlockSpec([], [:u_r, :u_i]))
        # slack node does not need to know the flowsum
        gen = generate_io_function(iob;
                                   f_states=[iob.u_r, iob.u_i],
                                   f_inputs=[],
                                   f_params=p_order,
                                   type=:ode)
        f = CallableBlockWrapper(gen)
    else
        error("The block should have outputs u_r and u_i. The inputs should be either i_r and i_i or empty (in case of a slack).")
    end

    vars = Symbol.(gen.states)
    ODEVertex(; f, dim=length(vars), sym=vars, mass_matrix=gen.massm)
end

# special function definition for empty bus
# function BusBar(; name=gensym(:Bus), verbose=false)
#     # busBar without any elements
#     @parameters t i_r(t) i_i(t) C ω0
#     @variables u_r(t) u_i(t)

#     bar = IOBlock([dt(u_r) ~ -ω0 * u_i - 1 / C * i_r,
#                    dt(u_i) ~ -ω0 * u_r - 1 / C * i_i],
#                   [i_r, i_i], [u_r, u_i];
#                   iv=t, name, warn=false)
# end

export BusBar
function BusBar(injectors...; name=:Bus, verbose=false, autopromote=false, EMT=true)
    for inj in injectors
        @assert BlockSpec([:u_r, :u_i], [:i_r, :i_i]; in_strict=false)(inj) "Injector $inj does not satisfy injector interface!"
    end

    # TODO assert that all iv ar equal, prob done in BlockSystems?
    @variables t
    ir, ii = Num[], Num[]
    for i in 1:length(injectors)
        irs = subscript(:i_r, i)
        iqs = subscript(:i_i, i)
        append!(ir, @parameters $irs(t))
        append!(ii, @parameters $iqs(t))
    end

    @parameters i_r(t) i_i(t)
    @variables u_r(t) u_i(t)
    dt = Differential(t)

    addall(list...) = isempty(list) ? 0 : (+)(list...)

    if EMT
        @parameters C ω0
        # @named bar = IOBlock([dt(u_r) ~  ω0*u_i - 1/C * (addall(ir...) + i_r),
        #                       dt(u_i) ~ -ω0*u_r - 1/C * (addall(ii...) + i_i)],
        @named bar = IOBlock([dt(u_r) ~ -ω0*u_i + 1/C * (addall(ir...) + i_r),
                              dt(u_i) ~  ω0*u_r + 1/C * (addall(ii...) + i_i)],
                             [i_r, i_i, ir..., ii...],
                             [u_r, u_i];
                             iv=t, warn=false)
    else
        @named bar = IOBlock([0 ~ (addall(ir...) + i_r),
                              0 ~ (addall(ii...) + i_i)],
                             [i_r, i_i, ir..., ii...],
                             [u_r, u_i];
                             iv=t, warn=false)
    end

    connections = Pair[]
    for (i, inj) in enumerate(injectors)
        push!(connections, inj.i_r => getproperty(bar, subscript(:i_r, i)))
        push!(connections, inj.i_i => getproperty(bar, subscript(:i_i, i)))
        push!(connections, bar.u_r => inj.u_r)
        push!(connections, bar.u_i => inj.u_i)
    end

    promotions = [bar.u_r => :u_r,
                  bar.u_i => :u_i,
                  bar.i_r => :i_r,
                  bar.i_i => :i_i]

    if EMT
        promotions = vcat(promotions, [bar.ω0=>:ω0, bar.C=>:C])
    end

    sys = IOSystem(connections, [bar, injectors...];
                   namespace_map=promotions, autopromote,
                   outputs=[bar.u_r, bar.u_i],
                   name)

    return connect_system(sys; verbose=verbose)
end

"""
Conventions:
(Anti-)symmetric lines:
  - inputs:   `u_src_r`, `u_src_i`, `u_dst_r`, `u_dst_i`
  - outputs:  `i_r`, `i_i`
  - current direction is defined from src to dst
  - dst node will receive `-i_r`, `-i_i`
Asymmetric lines:
  - inputs:   `u_src_r`, `u_src_i`, `u_dst_r`, `u_dst_i`
  - outputs:  `i_src_r`, `i_src_i`, `i_dst_r`, `i_dst_i`
  - not yet implemented. might by tricky with fidutial...
"""
function NetworkDynamics.ODEEdge(iob::IOBlock, p_order=[])
    symspec = BlockSpec([:u_src_r, :u_src_i, :u_dst_r, :u_dst_i],
                        [:i_r, :i_i])
    asymspec = BlockSpec([:u_src_r, :u_src_i, :u_dst_r, :u_dst_i],
                         [:i_src_r, :i_src_i, :i_dst_r, :i_dst_i])

    if fulfills(iob, symspec)
        return _odeedge_symmetric(iob::IOBlock, p_order)
    elseif fulfills(iob, asymspec)
        error("asymetric edges not yet implementet")
        return _odeedge_asymmetric(iob::IOBlock, p_order)
    else
        error("Unsupportet input-output relation in edge block.")
    end
end

function _odeedge_symmetric(iob, p_order)
    gen = generate_io_function(iob;
                               f_states=[iob.i_r, iob.i_i],
                               f_inputs=[iob.u_src_r, iob.u_src_i,
                                         iob.u_dst_r, iob.u_dst_i],
                               f_params=p_order,
                               warn=true,
                               type=:ode)
    N = length(gen.states)
    N != 2 && error("Does not work with internal states in the line right know!")

    f = function (du, u, src, dst, p, t)
        @views gen.f_ip(du[1:2], u[1:2], (src[1], src[2], dst[1], dst[2]), p, t)
        @views gen.f_ip(du[3:4], u[3:4], (dst[1], dst[2], src[1], src[2]), p, t)
    end
    names = String.(Symbol.(gen.states))
    sym = Symbol.(vcat(names .* "_dst", names .* "_src"))

    ODEEdge(; f, dim=4, coupling=:fiducial, mass_matrix=gen.massm, sym)
end

function NetworkDynamics.StaticEdge(iob::IOBlock, p_order=[])
    symspec = BlockSpec([:u_src_r, :u_src_i, :u_dst_r, :u_dst_i],
                        [:i_r, :i_i])
    asymspec = BlockSpec([:u_src_r, :u_src_i, :u_dst_r, :u_dst_i],
                         [:i_src_r, :i_src_i, :i_dst_r, :i_dst_i])

    if fulfills(iob, symspec)
        error("asymetric edges not yet implementet")
        return _staticedge_symmetric(iob::IOBlock, p_order)
    elseif fulfills(iob, asymspec)
        return _staticedge_asymmetric(iob::IOBlock, p_order)
    else
        error("Unsupportet input-output relation in edge block.")
    end
end

function _staticedge_asymmetric(iob::IOBlock, p_order)
    let gen = generate_io_function(iob;
                                   f_states=[iob.i_dst_r,
                                             iob.i_dst_i,
                                             iob.i_src_r,
                                             iob.i_src_i],
                                   f_inputs=[iob.u_src_r,
                                             iob.u_src_i,
                                             iob.u_dst_r,
                                             iob.u_dst_i],
                                   f_params=p_order,
                                   warn=true,
                                   type=:static)
        f = function (u, src, dst, p, t)
            gen.f_ip(u, (src[1], src[2], dst[1], dst[2]), p, t)
        end
        dim = length(gen.states)
        @assert dim == 4
        StaticEdge(; f, dim, coupling=:fiducial, sym=[:i_dst_r, :i_dst_i, :i_src_r, :i_src_i])
    end
end

function flowsum(edges)
    i_r, i_i = 0.0, 0.0
    for e in edges
        i_r += e[1]
        i_i += e[2]
    end
    return (i_r, i_i)
end


include("Components.jl")
include("inverter.jl")
include("models.jl")
include("solution_inspection.jl")
include("SerializeDict.jl")
include("import_network.jl")
include("interactive.jl")

end
