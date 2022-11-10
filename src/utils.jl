using Unitful

# define perunit as unit
@unit pu "p.u." PerUnit 1 false
Unitful.register(VirtualInertia)
function __init__()
    Unitful.register(VirtualInertia)
end

"""
    subscript(s, i, add...)

Append symbol or string `s` with a integer subscript.
"""
function subscript(s, i, add...)
    Symbol(s, _subscript(i), add...)
end

function _subscript(i::Int)
    dig = reverse(digits(i))
    String(map(i -> Char(0x02080 + i), dig))
end

# dq transformations
export Tdq, Tdqinv

Tdq(δ) = √(2/3)*[ cos(δ)  cos(δ-2pi/3)  cos(δ+2pi/3)
                 -sin(δ) -sin(δ-2pi/3) -sin(δ+2pi/3)]

Tdqinv(δ) = √(2/3)* [cos(δ)       -sin(δ)
                     cos(δ-2pi/3) -sin(δ-2pi/3)
                     cos(δ+2pi/3) -sin(δ+2pi/3)]

Tdq(t, a, b, c; ω=2π*50) = Tdq(ω*t)*[a, b, c]
function Tdqinv(t, d, q; ω=2π*50)
    N = length(t)
    a = similar(d)
    b = similar(d)
    c = similar(d)

    for i in 1:N
        a[i], b[i], c[i] = Tdqinv(ω*t[i])*[d[i], q[i]]
    end

    return a, b, c
end

export fixdiv
function fixdiv(A)
    for idx in eachindex(A)
        if A[idx] isa SymbolicUtils.Div
            A[idx] = eval(Meta.parse(repr(A[idx])))
        end
    end
    A = BlockSystems.narrow_type(A)
end

export u0guess
u0guess(nd::ODEFunction) = u0guess.(nd.syms)
function u0guess(s::Symbol)
    s = string(s)
    if occursin(r"^u_r", s)
        1.0
    elseif occursin(r"^A", s)
        1.0
    elseif occursin(r"^MfIf", s)
        sqrt(2/3) * 1/(2π*50)
    else
        0.0
    end
end

treesyms(x) = treesyms(length(x))
treesyms(i::Integer) = Iterators.flatten((Iterators.repeated('├', i-1), '└'))

export PRecord, record!
struct PRecord{PT}
    t::Vector{Float64}
    p::Vector{PT}
end

foo(x) ="bar"
PRecord(prob::ODEProblem) = PRecord([Float64(prob.tspan[1])], [prob.p])

function record!(pr::PRecord, integrator)
    push!(pr.t, integrator.t)
    push!(pr.p, deepcopy(integrator.p))
end

function (pr::PRecord)(t::Number; direction=:right)
    fun = direction==:right ? tr->tr>t : tr->tr≥t
    idx = findfirst(fun, pr.t)
    if isnothing(idx)
        pr.p[end]
    elseif idx == 1
        pr.p[begin]
    else
        pr.p[idx - 1]
    end
end

function (pr::PRecord{PT})(ts) where {PT}
    p = Vector{PT}(undef, length(ts))
    lastt = -Inf
    for (i, t) in enumerate(ts)
        if t == lastt
            p[i] = pr(t; direction=:right)
        else
            p[i] = pr(t; direction=:left)
        end
        lastt = t
    end
    p
end
