export Components

module Components

using BlockSystems

export LowPassFilter, DroopControl, VoltageSource, Power, PowerConstraint, InversePowerConstraint, ImpedanceConstraint


"""
    LowPassFilter(;name, renamings...)

Returns a low pass filter. The name of the system and the names of the vars
can be changed with keyword arguments `name=:myname, τ=:mytau, …`.

    out'(t) = 1/τ (in(t) - out(t))

               +-----+
    input(t) --|  τ  |-- output(t)
               +-----+

    IOBlock :##LPF# with 1 eqs
    ├ inputs:  input(t)
    ├ outputs: output(t)
    ├ istates: (empty)
    └ iparams: τ
"""
function LowPassFilter(;name=gensym(:lpf), renamings...)
    @parameters t τ
    @parameters input(t)
    @variables output(t)
    D = Differential(t)

    block = IOBlock([D(output) ~ 1/τ * (- output + input)],
                [input], [output]; name)
    return isempty(renamings) ? block : rename_vars(block; renamings...)
end

"""
    Droop(;name, renamings...)

Returns a Droop. The name of the system and the names of the vars
can be changed with keyword arguments `name=:myname, K=:myK, …`.

    u = - K*(x - x_ref) + u_ref

           +-----------------+
    x(t) --| K, x_ref, u_ref |-- u(t)
           +-----------------+

    IOBlock :##droop# with 1 eqs
    ├ inputs:  x(t)
    ├ outputs: u(t)
    ├ istates: (empty)
    └ iparams: K, x_ref, u_ref
"""
function Droop(;name=gensym(:droop), renamings...)
    @parameters t K x_ref u_ref
    @parameters x(t)
    @variables u(t)
    D = Differential(t)

    block = IOBlock([u ~ - K * (x - x_ref) + u_ref], # output is the droop voltage v
                    [x], [u]; name)

    return isempty(renamings) ? block : rename_vars(block; renamings...)
end

"""
    Power(;name, renamings...)

Returns a Block which calculates the active and reactive power for a given complex input.

    P = uᵣ iᵣ + uᵢ iᵢ
    Q = uᵢ iᵣ - uᵣ iᵢ

             +-----+
    u_r(t) --|     |-- P(t)
    u_i(t) --|     |
    i_r(t) --|     |
    i_i(t) --|     |-- Q(t)
             +-----+

    IOBlock :##power# with 2 eqs
    ├ inputs:  u_i(t), u_r(t), i_i(t), i_r(t)
    ├ outputs: P(t), Q(t)
    ├ istates: (empty)
    └ iparams: (empty)
"""
function Power(;name=gensym(:power), renamings...)
    @parameters t
    @parameters u_i(t) u_r(t) i_i(t) i_r(t)
    @variables P(t) Q(t)

    block = IOBlock([P ~ u_r*i_r + u_i*i_i,
                     Q ~ u_i*i_r - u_r*i_i],
                    [u_i, u_r, i_i, i_r], [P, Q]; name)

    return isempty(renamings) ? block : rename_vars(block; renamings...)
end

"""
    Cart2Polar(;name=:c2p, renamings...)

(X, Y) ↦ (mag, arg) transformation
"""
function Cart2Polar(;name=:c2p, renamings...)
    @variables t arg(t) mag(t)
    @parameters x(t) y(t)
    block = IOBlock([mag ~ √(x^2 + y^2),
                     arg ~ atan(y, x)],
                    [x, y], [mag, arg]; name)

    rename_vars(block; renamings...)
end
"""
    Polar2Cart(;name=:p2c, renamings...)

(mag, arg) ↦ (X, Y) transformation
"""
function Polar2Cart(;name=:p2c, renamings...)
    @variables t x(t) y(t)
    @parameters arg(t) mag(t)
    block = IOBlock([x ~ mag * cos(arg),
                     y ~ mag * sin(arg)],
                    [mag, arg], [x, y]; name)

    rename_vars(block; renamings...)
end

"""
    VRefGen(; name=:Vrefgen, renamings...)

Create dq-reference from (ω, V) reference.
"""
function VRefGen(; name=:Vrefgen, renamings...)
    @variables t δ(t) u_r_ref(t) u_i_ref(t)
    @parameters V(t) ω(t)
    dt = Differential(t)
    block = IOBlock([dt(δ) ~ ω,
                     u_r_ref ~ V*cos(δ),
                     u_i_ref ~ V*sin(δ)],
                    [V, ω],
                    [u_r_ref, u_i_ref];
                    name)
    rename_vars(block; renamings...)
end

"""
    UIMeas(; name=:uimeas, renamings...)

Creates a very simple block which is mainly there for renaiming.
"""
function UIMeas(; name=:ui_meas, renamings...)
    @variables t u_r_meas(t) u_i_meas(t) i_r_meas(t) i_i_meas(t)
    @parameters u_r(t) u_i(t) i_i(t) i_r(t)
    block = IOBlock([u_r_meas ~ u_r,
                     u_i_meas ~ u_i,
                     i_r_meas ~ -i_r,
                     i_i_meas ~ -i_i],
                    [u_r, u_i, i_r, i_i],
                    [u_r_meas, u_i_meas, i_r_meas, i_i_meas];
                    name)

    rename_vars(block; renamings...)
end

function ReducedPLL(; name=:pll, renamings...)
    @variables t δ_pll(t) ω_pll(t) μ_pll(t) u_qpll(t) u_q(t)
    @parameters u_r(t) u_i(t) Kp Ki ω_lp
    dt = Differential(t)
    block = IOBlock([u_q ~ 1/sqrt(u_r^2+u_i^2)*(u_r*sin(-δ_pll) + u_i*cos(-δ_pll)),
                     dt(u_qpll) ~ ω_lp*(u_q - u_qpll),
                     dt(μ_pll) ~ u_qpll,
                     ω_pll ~ Kp*u_qpll + Ki*μ_pll,
                     dt(δ_pll) ~ ω_pll
                     ],
                    [u_r, u_i],
                    [δ_pll, ω_pll];
                    name)
    block = substitute_algebraic_states(block)
    rename_vars(block; renamings...)
end

function KauraPLL(; name=:pll, renamings...)
    @variables t δ_pll(t) ω_pll(t) ϵ_pll(t) u_d_pll(t) u_q_pll(t) u_d_out(t) u_q_out(t)
    @parameters u_r(t) u_i(t) Kp Ki ω_lp
    dt = Differential(t)
    block = IOBlock([u_d_out ~ 1/sqrt(u_r^2+u_i^2)*(u_r*cos(-δ_pll) - u_i*sin(-δ_pll)),
                     u_q_out ~ 1/sqrt(u_r^2+u_i^2)*(u_r*sin(-δ_pll) + u_i*cos(-δ_pll)),
                     dt(u_d_pll) ~ ω_lp*(u_d_out - u_d_pll),
                     dt(u_q_pll) ~ ω_lp*(u_q_out - u_q_pll),
                     dt(ϵ_pll) ~ atan(u_q_pll, u_d_pll),
                     ω_pll ~ Kp*atan(u_q_pll, u_d_pll) + Ki* ϵ_pll,
                     dt(δ_pll) ~ ω_pll
                     ],
                    [u_r, u_i],
                    [δ_pll, ω_pll];
                    name)
    block = substitute_algebraic_states(block)
    rename_vars(block; renamings...)
end

function PT1CurrentSource_old(; name=:CS, renamings...)
    @variables t i_r(t) i_i(t)
    @parameters P_ref(t) τ u_r(t) u_i(t)
    dt = Differential(t)
    block = IOBlock([dt(i_r) ~ 1/τ * (P_ref * u_r/(u_r^2 + u_i^2) - i_r),
                     dt(i_i) ~ 1/τ * (P_ref * u_i/(u_r^2 + u_i^2) - i_i)],
                    [u_r, u_i, P_ref], [i_r, i_i],
                    name=:CS)
    replace_vars(block; renamings...)
end

function PT1CurrentSource(; name=:CS, renamings...)
    lpf_r = Components.LowPassFilter(; input=:i_ref_r, output=:i_r)
    lpf_i = Components.LowPassFilter(; input=:i_ref_i, output=:i_i)
    pcs   = Components.PerfectCurrentSource(; i_r=:i_ref_r, i_i=:i_ref_i)

    sys = IOSystem(:autocon, [lpf_r, lpf_i, pcs]; globalp=[:τ], outputs=:remaining, name)
    block = connect_system(sys)
    replace_vars(block; renamings...)
end

function PerfectCurrentSource(; name=:CS, renamings...)
    @variables t i_r(t) i_i(t)
    @parameters u_r(t) u_i(t) P Q
    block = IOBlock([i_r ~ (P*u_r + Q*u_i)/(u_r^2 + u_i^2),
                     i_i ~ (P*u_i - Q*u_r)/(u_r^2 + u_i^2)],
                    [u_r, u_i], [i_r, i_i],
                    name=:CS)
    replace_vars(block; renamings...)
end

end # module
