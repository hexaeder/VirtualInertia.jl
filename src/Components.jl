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

end # module
