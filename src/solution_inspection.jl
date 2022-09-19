using LinearAlgebra: norm
using SciMLBase
export blockstates, getstate, timeseries


blockstates(sol, idx) = blockstates(sol.prob.f, idx)
function blockstates(nd::ODEFunction, idx)
    wrapper = _getwrapper(nd,idx)
    vcat(wrapper.states, wrapper.rem_states)
end

function _getwrapper(nd, idx)
    ndobj = nd.f
    _group = findfirst(group -> idx ∈ group, ndobj.unique_v_indices)
    ndobj.unique_vertices![_group].f
end

getstate(sol, t::Number, idx, state) = getstate(sol, t, nothing, idx, state)
function getstate(sol, t::Number, p, idx, state)
    nd = sol.prob.f
    x = sol(t)
    gd = nd(x, p, t, GetGD)
    vstate = collect(get_vertex(gd, idx))
    wrapper = _getwrapper(nd, idx)

    if state ∈ wrapper.states
        stateidx = findfirst(s->s==state, wrapper.states)
        return vstate[stateidx]
    elseif state ∈ wrapper.rem_states
        stateidx = findfirst(s->s==state, wrapper.rem_states)
        input = flowsum(get_dst_edges(gd, idx))
        pblock = p isa Tuple ? p[1][idx] : nothing
        return wrapper.g_oop(vstate, input, pblock, t)[stateidx]
    elseif state==:Vmag
        return norm(vstate[1:2])
    elseif state==:Varg
        return atan(vstate[2], vstate[1])
    elseif state==:imag
        input = flowsum(get_dst_edges(gd, idx))
        return norm(input[1:2])
    elseif state==:iarg
        input = flowsum(get_dst_edges(gd, idx))
        return atan(input[2], input[1])
    elseif state==:Vmag_ref
        u_r_ref = getstate(sol, t, p, idx, :u_r_ref)
        u_i_ref = getstate(sol, t, p, idx, :u_i_ref)
        return norm([u_r_ref, u_i_ref])
    elseif state==:Varg_ref
        u_r_ref = getstate(sol, t, p, idx, :u_r_ref)
        u_i_ref = getstate(sol, t, p, idx, :u_i_ref)
        return atan(u_i_ref, u_r_ref)
    elseif state==:ωmeas
        u_r, u_i = vstate[1:2]
        dx = sol(t, Val{1})
        gd = nd(dx, p, t, GetGD)
        dvstate = collect(get_vertex(gd, idx))
        u_r_dot, u_i_dot = dvstate[1:2]
        return -(u_i*u_r_dot - u_r*u_i_dot)/(u_i^2 + u_r^2)
    elseif state==:Va
        return (Tdqinv(2π*50*t)*[vstate[1], vstate[2]])[1]
    elseif state==:Vb
        return (Tdqinv(2π*50*t)*[vstate[1], vstate[2]])[2]
    elseif state==:Vc
        return (Tdqinv(2π*50*t)*[vstate[1], vstate[2]])[3]
    elseif state==:ia
        i_r, i_i = flowsum(get_dst_edges(gd, idx))
        return (Tdqinv(2π*50*t)*[i_r, i_i])[1]
    elseif state==:ib
        i_r, i_i = flowsum(get_dst_edges(gd, idx))
        return (Tdqinv(2π*50*t)*[i_r, i_i])[2]
    elseif state==:ic
        i_r, i_i = flowsum(get_dst_edges(gd, idx))
        return (Tdqinv(2π*50*t)*[i_r, i_i])[3]
    elseif state==:Smeas
        u_r, u_i = vstate[1:2]
        i_r, i_i = flowsum(get_dst_edges(gd, idx))
        return (u_r + im*u_i)*(-i_r + im*i_i)
    elseif state==:Pmeas
        return real(getstate(sol, t, p, idx, :Smeas))
    elseif state==:Qmeas
        return imag(getstate(sol, t, p, idx, :Smeas))
    else
        error("Don't know state $state.")
    end
end

timeseries(sol, ts, idx::Int, state) = timeseries(sol, ts, nothing, idx, state)
timeseries(sol, ts, p, idx::Int, state) = (collect(ts), [getstate(sol, t, p, idx, state) for t in ts])

TS_DTMAX::Float64 = 0.01
set_ts_dtmax(dt) = global TS_DTMAX=dt
timeseries(sol, idx::Int, state; kwargs...) = timeseries(sol, nothing, idx, state; kwargs...)
function timeseries(sol, p, idx::Int, state; dtmax=TS_DTMAX)
    if p==nothing && sol.prob.p !==nothing
        p = sol.prob.p
        if p !== SciMLBase.NullParameters()
            @warn "I am using the p from the problem $p to recover states. Be carefull, this might be wrong afer callbacks."
        end
    end

    ts = [sol.t[begin]]
    for i in eachindex(sol.t)[begin+1:end]
        while !isnothing(dtmax) && ts[end] + dtmax < sol.t[i]
            push!(ts, ts[end]+dtmax)
        end
        push!(ts, sol.t[i])
    end
    timeseries(sol, ts, p, idx, state)
end
