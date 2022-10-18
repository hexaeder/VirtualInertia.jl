using LinearAlgebra: norm
using SciMLBase
export blockstates, getstate, timeseries, meanseries
export NADIR, ROCOF, needed_storage

const TS_DTMAX = Ref(0.01)

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

@nospecialize
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
    elseif state==:rocof
        h = 0.005
        t1 = t-h < sol.t[begin] ? t : t-h
        t2 = t+h < sol.t[begin] ? t : t+h
        ω1 = getstate(sol, t1, p, idx, :ωmeas)
        ω2 = getstate(sol, t2, p, idx, :ωmeas)
        return (ω2-ω1)/(t2-t1)
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

_timeseries(sol, ts, idx::Int, state) = _timeseries(sol, ts, nothing, idx, state)
_timeseries(sol, ts, p, idx::Int, state) = (collect(ts), [getstate(sol, t, p, idx, state) for t in ts])

set_ts_dtmax(dt) = VirtualInertia.TS_DTMAX[] = dt
timeseries(sol, idx::Int, state; kwargs...) = timeseries(sol, nothing, idx, state; kwargs...)
function timeseries(sol, p, idx::Int, state; dtmax=TS_DTMAX[])
    if p==nothing && sol.prob.p !==nothing
        p = sol.prob.p
        if p !== SciMLBase.NullParameters()
            # @warn "I am using the p from the problem $p to recover states. Be carefull, this might be wrong afer callbacks."
        end
    end

    ts = Float64[]
    for i in eachindex(sol.t)
        if isempty(ts) || ts[end] !== sol.t[i]
            push!(ts, sol.t[i])
        end

        if i<lastindex(sol.t)
            tdiff = sol.t[i+1] - sol.t[i]
            if tdiff > dtmax
                nsteps = tdiff ÷ dtmax + 1
                dt = tdiff/nsteps
                for k in 1:nsteps-1
                    push!(ts, ts[end] + dt)
                end
            end
        end
    end
    _timeseries(sol, ts, p, idx, state)
end

function meanseries(sol, idxs, state; dtmax=TS_DTMAX[])
    ts, x = timeseries(sol, idxs[begin], state; dtmax)
    for i in idxs[begin+1:end]
        _, xi = _timeseries(sol, ts, i, state)
        x += xi
    end
    x = x/length(idxs)
    return ts, x
end

function NADIR(sol, idx::Int)
    _, ω = timeseries(sol, idx, :ωmeas, dtmax=0.01)
    return argmax(abs, ω)
end
function NADIR(sol, idxs)
    d = Dict{Any, Float64}(idxs .=> NADIR.(Ref(sol), idxs))
    _, ωmean = meanseries(sol, idxs, :ωmeas; dtmax=0.01)
    d[:mean] = argmax(abs, ωmean)
    d
end

function ROCOF(sol, idx::Int)
    _, rocof = timeseries(sol, idx, :rocof, dtmax=0.01)
    return argmax(abs, rocof)
end

function ROCOF(sol, idxs)
    d = Dict{Any, Float64}(idxs .=> ROCOF.(Ref(sol), idxs))
    _, rocofmean = meanseries(sol, idxs, :rocof; dtmax=0.01)
    d[:mean] = argmax(abs, rocofmean)
    d
end

needed_storage(sol, idxs) = needed_storage.(Ref(sol), idxs)
function needed_storage(sol, idx::Int)
    node_p = sol.prob.p[1][idx]
    psyms = VirtualInertia._getwrapper(sol.prob.f, idx).params
    pidx = findfirst(isequal(:P_ref), psyms)

    t, P = timeseries(sol, idx, :Pmeas)
    Pover = P .- node_p[pidx]
    sum(diff(t) .* Pover[1:end-1])
end

function powerloss(sol)
    t, P = timeseries(sol, 1, :Pmeas)
    tidx = findfirst(x->x ≥ 0.1, t)
    t = t[tidx:end]
    P = P[tidx:end]
    Punder = 0.1 .- P
    sum(diff(t) .* Punder[1:end-1])
end

export plotsym, plotsym!
import Plots
plotsym(sol::ODESolution, sym, nodes=1:nv(sol.prob.f.f.graph); mean=false) = plotsym!(Plots.plot(), sol, sym, nodes; mean)
plotsym!(sol::ODESolution, sym, nodes=1:nv(sol.prob.f.f.graph); mean=false) = plotsym!(Plots.current(), sol, sym, nodes; mean)
function plotsym!(p::Plots.Plot, sol::ODESolution, sym, nodes=1:nv(sol.prob.f.f.graph); mean=false)
    if mean
        Plots.plot!(p, meanseries(sol, nodes, sym); label=string(sym)*" mean")
    else
        for i in nodes
            try
                Plots.plot!(p, timeseries(sol, i, sym); label=string(sym)*string(i))
            catch
                @warn "Could not create timeseries for $sym on node $i. Skipped!"
            end
        end
    end
    p
end

@specialize
