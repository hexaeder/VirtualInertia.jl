using LinearAlgebra: norm
using SciMLBase
import Makie
import Makie.GeometryBasics: Point2
import RecipesBase
export blockstates, getstate, timeseries, meanseries
export NADIR, ROCOF, needed_storage

const TS_DTMAX = Ref(0.01)
const MEAS_STATES = [:_ω, :_rocof, :_P, :_Q]
const MEAS_STATES_COMPLEX = [:_u, :_i, :_S]
const COMPLEX_ENDINGS_CART = ["_r", "_i"]
const COMPLEX_ENDINGS_POL = ["_mag", "_arg"]

blockstates(sol::ODESolution, idx; print=true) = blockstates(sol.prob.f, idx; print)
function blockstates(nd::ODEFunction, idx; print=true, polar=true)
    wrp = _getwrapper(nd,idx)

    states = _group_states(wrp.states)
    remstates = _group_states(wrp.rem_states)
    measstates = [MEAS_STATES, MEAS_STATES_COMPLEX]

    if print
        println("Actual States")
        sstates = _stringify(states...)
        for (s, c) in zip(sstates, treesyms(sstates))
            println(" $c $s")
        end
        println("Removed States")
        sstates = _stringify(remstates...)
        for (s, c) in zip(sstates, treesyms(sstates))
            println(" $c $s")
        end
        println("Measured States (available for all nodes)")
        sstates = _stringify(measstates...)
        for (s, c) in zip(sstates, treesyms(sstates))
            println(" $c $s")
        end
        println("All complex states can be used with _r, _i, _mag, _arg, _a, _b, _c suffix")
        return nothing
    end

    cmplx_endings = if polar
        COMPLEX_ENDINGS_POL
    else
        COMPLEX_ENDINGS_CART
    end
    for s in (states, remstates, measstates)
        s[2] = Symbol.(sort!(mapreduce(t->t[1]*t[2], vcat, Iterators.product(string.(s[2]), cmplx_endings))))
        append!(s[2], s[1])
    end
    return (;states=states[2], rem_states=remstates[2], meas_states=measstates[2])
end

function _common_states(sol, idxs)
    idxs = collect(idxs)
    nt = blockstates(sol, idxs[1]; print=false)
    states = Set{Symbol}(nt.states)
    rem = Set{Symbol}(nt.rem_states)
    meas = nt.meas_states
    for i in idxs[2:end]
        nt = blockstates(sol, i; print=false)
        intersect!(states, nt.states)
        intersect!(rem, nt.rem_states)
    end
    return (; states=sort!(collect(states)), rem_states=sort!(collect(rem)), meas_states=meas)
end

_getwrapper(sol::ODESolution, idx) = _getwrapper(sol.prob.f, idx)
function _getwrapper(nd, idx)
    ndobj = nd.f
    _group = findfirst(group -> idx ∈ group, ndobj.unique_v_indices)
    ndobj.unique_vertices![_group].f
end

"""
Find complex states with _r and _i suffix
"""
function _group_states(states)
    states = string.(states)
    complex_pairs = String[]
    for state in states
        if occursin(r"_r$", state)
            m = match(r"^(.*)_r$", state)
            name = m[1]
            if !isnothing(findfirst(isequal(name*"_i"), states))
                push!(complex_pairs, name)
            end
        end
    end
    real = Symbol.(sort!(filter(s -> s ∉ complex_pairs .* "_r" ∪ complex_pairs .* "_i", states)))
    complex = Symbol.(complex_pairs)
    [real, complex]
end

function _stringify(real, complex)
    complex = replace.(string.(complex), r"^(.*)$" => s"Complex: \1")
    append!(complex, string.(real))
end



@nospecialize
function getstate(sol, t::Number, p, idxs, state)
   [getstate(sol, t, p, idx, state; err=false) for idx in idxs]
end

function getstate(sol, t::Number, p, idx::Integer, state; err=true)
    if p isa PRecord
        p = p(t)
    end
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
    # measured states below
    elseif state==:_u_r
        return vstate[1]
    elseif state==:_u_i
        return vstate[2]
    elseif state==:_i_r
        flowsum(get_dst_edges(gd, idx))[1]
    elseif state==:_i_i
        flowsum(get_dst_edges(gd, idx))[2]
    elseif state==:_i_mag
        input = flowsum(get_dst_edges(gd, idx))
        return norm(input[1:2])
    elseif state==:_i_arg
        input = flowsum(get_dst_edges(gd, idx))
        return atan(input[2], input[1])
    elseif occursin(r"_mag" , string(state))
        name = match(r"^(.*)_mag$", string(state))[1]
        _r = getstate(sol, t, p, idx, Symbol(name*"_r"); err)
        _i = getstate(sol, t, p, idx, Symbol(name*"_i"); err)
        return norm((_r, _i))
    elseif occursin(r"_arg" , string(state))
        name = match(r"^(.*)_arg$", string(state))[1]
        _r = getstate(sol, t, p, idx, Symbol(name*"_r"); err)
        _i = getstate(sol, t, p, idx, Symbol(name*"_i"); err)
        return atan(_i, _r)
    elseif state==:_ω
        u_r, u_i = vstate[1:2]
        dx = sol(t, Val{1})
        gd = nd(dx, p, t, GetGD)
        dvstate = collect(get_vertex(gd, idx))
        u_dot_r, u_dot_i = dvstate[1:2]
        return -(u_i*u_dot_r - u_r*u_dot_i)/(u_i^2 + u_r^2)
    elseif state==:_rocof
        h = 0.001
        # t1 = t-h < sol.t[begin] ? t : t-h
        t1 = t
        t2 = t+h < sol.t[begin] ? t : t+h
        ω1 = getstate(sol, t1, p, idx, :_ω; err)
        ω2 = getstate(sol, t2, p, idx, :_ω; err)
        return (ω2-ω1)/(t2-t1)
    elseif occursin(r"_[abc]" , string(state))
        m = match(r"^(.*)_(.)$", string(state))
        name = m[1]
        phase = m[2]
        phase = Dict("a"=>1, "b"=>2, "c"=>3)[phase]
        _r = getstate(sol, t, p, idx, Symbol(name*"_r"); err)
        _i = getstate(sol, t, p, idx, Symbol(name*"_i"); err)
        return (Tdqinv(2π*50*t)*[vstate[1], vstate[2]])[phase]
    elseif state==:_S
        u_r, u_i = vstate[1:2]
        i_r, i_i = flowsum(get_dst_edges(gd, idx))
        return (u_r + im*u_i)*(-i_r + im*i_i)
    elseif state==:_P || state==:_S_r
        return real(getstate(sol, t, p, idx, :_S; err))
    elseif state==:_Q || state==:_S_i
        return imag(getstate(sol, t, p, idx, :_S; err))
    else
        if err
            error("Don't know state $state. Call `blocksates(sol, idx)` to get list of states.")
        else
            NaN
        end
    end
end

struct TimeSeries
    t::Vector{Float32}
    x::Vector{Float32}
    name::String
end

Base.getindex(ts::TimeSeries, i) = ts.x[i]
Base.length(ts::TimeSeries) = length(ts.x)
Base.firstindex(ts::TimeSeries) = firstindex(ts.x)
Base.lastindex(ts::TimeSeries) = lastindex(ts.x)

function _blockname(sol, idx)
    try
        return " ("*_getwrapper(sol, idx).name*")"
    catch
    end
    return ""
end

_timeseries(sol, p::PRecord, ts, idx::Int, state) = TimeSeries(collect(ts), [Float32(getstate(sol, t, p(t), idx, state)) for t in ts], "$state @ $idx"*_blockname(sol,idx))
_timeseries(sol, p, ts, idx::Int, state) = TimeSeries(collect(ts), [Float32(getstate(sol, t, p, idx, state)) for t in ts], "$state @ $idx"*_blockname(sol,idx))

set_ts_dtmax(dt) = VirtualInertia.TS_DTMAX[] = dt
timeseries(sol, idx::Int, state; kwargs...) = timeseries(sol, nothing, idx, state; kwargs...)
function timeseries(sol, p, idx::Int, state; dtmax=TS_DTMAX[])
    if p==nothing && sol.prob.p !==nothing
        p = sol.prob.p
        if p !== SciMLBase.NullParameters()
            @warn "I am using the p from the problem to recover states. Be carefull, this might be wrong afer callbacks."
        end
    end

    ts = Float32[]
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
    _timeseries(sol, p, ts, idx, state)
end

function meanseries(sol, idxs, state; dtmax=TS_DTMAX[])
    ts, x = timeseries(sol, idxs[begin], state; dtmax)
    for i in idxs[begin+1:end]
        _, xi = _timeseries(sol, ts, i, state)
        x += xi
    end
    x = x/length(idxs)
    return Timeseries(ts, x, "$state mean")
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
                Plots.plot!(p, timeseries(sol, i, sym))
            catch
                @warn "Could not create timeseries for $sym on node $i. Skipped!"
            end
        end
    end
    p
end

@specialize

function Makie.convert_arguments(P::Makie.PointBased, ts::TimeSeries)
    return (Point2f.(zip(ts.t, ts.x)),)
end
Makie.plottype(::TimeSeries) = Makie.Lines

RecipesBase.@recipe function f(ts::TimeSeries)
    seriestype  :=  :path
    label   --> ts.name
    (ts.t, ts.x)
end
