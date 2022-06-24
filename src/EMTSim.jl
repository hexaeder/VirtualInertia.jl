module EMTSim

using NetworkDynamics
using BlockSystems
using ModelingToolkit

include("utils.jl")

function NetworkDynamics.ODEVertex(iob::IOBlock, p_order=[])
    @assert length(p_order) == length(iob.iparams) "Provide order of all iparams!"

    if fulfills(iob, BlockSpec([:i_d, :i_q], [:u_d, :u_q]))
        # normal node
        gen = generate_io_function(iob;
                                   f_states=[iob.u_d, iob.u_q],
                                   f_inputs=[iob.i_d, iob.i_q],
                                   f_params=p_order,
                                   warn=true,
                                   type=:ode)

        f = (du, u, edges, p, t) -> gen.f_ip(du, u, flowsum(edges), p, t)
    elseif fulfills(iob, BlockSpec([], [:u_d, :u_q]))
        # slack node does not need to know the flowsum
        gen = generate_io_function(iob;
                                   f_states=[iob.u_d, iob.u_q],
                                   f_inputs=[],
                                   f_params=p_order,
                                   warn=true,
                                   type=:ode)
        f = (du, u, edges, p, t) -> gen.f_ip(du, u, nothing, p, t)
    else
        error("The block should have outputs u_d and u_q. The inputs should be either i_d and i_q or empty (in case of a slack).")
    end

    vars = Symbol.(gen.states)
    ODEVertex(; f, dim=length(vars), sym=vars, mass_matrix=gen.massm)
end

# special function definition for empty bus
# function BusBar(; name=gensym(:Bus), verbose=false)
#     # busBar without any elements
#     @parameters t i_d(t) i_q(t) C ω0
#     @variables u_d(t) u_q(t)

#     bar = IOBlock([dt(u_d) ~ -ω0 * u_q - 1 / C * i_d,
#                    dt(u_q) ~ -ω0 * u_d - 1 / C * i_q],
#                   [i_r, i_i], [u_r, u_i];
#                   iv=t, name, warn=false)
# end

export BusBar
function BusBar(injectors...; name=gensym(:Bus), verbose=false)
    for inj in injectors
        @assert BlockSpec([:u_d, :u_q], [:i_d, :i_q])(inj) "Injector $inj does not satisfy injector interface!"
    end

    # TODO assert that all iv ar equal, prob done in BlockSystems?
    @variables t
    id, iq = Num[], Num[]
    for i in 1:length(injectors)
        ids = subscript(:i_d, i)
        iqs = subscript(:i_q, i)
        append!(id, @parameters $ids(t))
        append!(iq, @parameters $iqs(t))
    end

    @parameters i_d(t) i_q(t) C ω0
    @variables u_d(t) u_q(t)
    dt = Differential(t)

    addall(list...) = isempty(list) ? 0 : (+)(list...)

    @named bar = IOBlock([dt(u_d) ~  ω0*u_q + 1/C * (addall(id...) + i_d),
                          dt(u_q) ~ -ω0*u_d + 1/C * (addall(iq...) + i_q)],
                         [i_d, i_q, id..., iq...],
                         [u_d, u_q];
                         iv=t, warn=false)

    connections = Pair[]
    for (i, inj) in enumerate(injectors)
        push!(connections, inj.i_d => getproperty(bar, subscript(:i_d, i)))
        push!(connections, inj.i_q => getproperty(bar, subscript(:i_q, i)))
        push!(connections, bar.u_d => inj.u_d)
        push!(connections, bar.u_q => inj.u_q)
    end

    promotions = [bar.u_d => :u_d,
                  bar.u_q => :u_q,
                  bar.i_d => :i_d,
                  bar.i_q => :i_q,
                  bar.ω0  => :ω0,
                  bar.C   => :C]

    sys = IOSystem(connections, [bar, injectors...];
                   namespace_map=promotions, autopromote=false,
                   outputs=[bar.u_d, bar.u_q],
                   name)

    return connect_system(sys; verbose=verbose)
end

"""
Conventions:
(Anti-)symmetric lines:
  - inputs:   `u_d_src`, `u_q_src`, `u_d_dst`, `u_q_dst`
  - outputs:  `i_d`, `i_q`
  - current direction is defined from src to dst
  - dst node will receive `-i_d`, `-i_q`
Asymmetric lines:
  - inputs:   `u_d_src`, `u_q_src`, `u_d_dst`, `u_q_dst`
  - outputs:  `i_d_src`, `i_q_src`, `i_d_dst`, `i_q_dst`
  - not yet implemented. might by tricky with fidutial...
"""
function NetworkDynamics.ODEEdge(iob::IOBlock, p_order=[])
    symspec = BlockSpec([:u_d_src, :u_q_src, :u_d_dst, :u_q_dst],
                        [:i_d, :i_q])
    asymspec = BlockSpec([:u_d_src, :u_q_src, :u_d_dst, :u_q_dst],
                         [:i_d_src, :i_q_src, :i_d_dst, :i_q_dst])

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
                               f_states=[iob.i_d, iob.i_q],
                               f_inputs=[iob.u_d_src, iob.u_q_src,
                                         iob.u_d_dst, iob.u_q_dst],
                               f_params=p_order,
                               warn=true,
                               type=:ode)
    N = length(gen.states)
    N != 2 && error("Does not work with internal states in the line right know!")

    f = function (du, u, src, dst, p, t)
        # TODO: allocation free stack
        gen.f_ip(du, u, vcat(src, dst), p, t)
        du[(N + 1):(2N)] .= -1.0 .* du[1:N]
        # # du[(N + 1):(2N)] .= du[N + 1]
        # # du[1:N] .*= -1.0
        # du .= -1.0 .* du
        # println(du)
    end
    names = String.(Symbol.(gen.states))
    sym = Symbol.(vcat(names .* "_dst", names .* "_src"))

    ODEEdge(; f, dim=4, coupling=:fiducial, mass_matrix=gen.massm, sym)
end

function flowsum(edges)
    i_d, i_q = 0.0, 0.0
    for e in edges
        i_d += e[1]
        i_q += e[2]
    end
    return (i_d, i_q)
end

end
