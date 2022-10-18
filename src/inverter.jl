export Inverter, VSwithLoad, DroopControl, Synchronverter, PerfectSource, PT1Source

function Inverter(inner, outer; name=Symbol(string(outer.name)*"_"*string(inner.name)))
    @assert BlockSpec([:i_i, :i_r, :u_r_ref, :u_i_ref], [:u_r, :u_i])(inner) "Inner ctrl loop does not meet expectation."
    @assert BlockSpec([], [:u_r_ref, :u_i_ref]; in_strict=false)(outer)  "Outer ctrl loop does not meet expectation."

    outerinputs = ModelingToolkit.getname.(outer.inputs)
    :u_r ∈ outerinputs && @error "Outer shoud use :u_r_meas instead of :u_r"
    :u_i ∈ outerinputs && @error "Outer shoud use :u_i_meas instead of :u_i"
    :i_r ∈ outerinputs && @error "Outer shoud use :i_r_meas instead of :i_r"
    :i_i ∈ outerinputs && @error "Outer shoud use :i_i_meas instead of :i_i"

    sys = IOSystem(:autocon, [outer, inner];
                   outputs=:remaining,
                   name)
    closed = connect_system(sys)

    @assert BlockSpec([:i_i, :i_r], [:u_r, :u_i])(closed) "Closed loop does not match expectation! $closed"

    return closed
end

function VSwithLoad(vs; name=Symbol(string(vs.name)*"_w_load"), params...)
    @assert BlockSpec([:i_i, :i_r], [:u_r, :u_i])(vs) "Does not look like a voltage source: $closed"

    pcs = Components.PerfectCurrentSource(; name=:load, P=:P_load, Q=:Q_load, i_i=:i_load_i, i_r=:i_load_r)

    @variables t i_int_r(t) i_int_i(t)
    @parameters i_r(t) i_i(t) i_load_r(t) i_load_i(t)
    kirchhoff = IOBlock([i_int_r ~ i_r + i_load_r,
                         i_int_i ~ i_i + i_load_i],
                        [i_r, i_i, i_load_r, i_load_i], [i_int_r, i_int_i])

    vs = replace_vars(vs; i_r=:i_int_r, i_i=:i_int_i)

    sys = IOSystem(:autocon, [vs, pcs, kirchhoff]; outputs=[vs.u_r, vs.u_i], name)
    con = connect_system(sys)
end

####
#### Outer Control Loops
####
"""
    DroopControl(; params...)

Return block for droop control outer ocntrol.
"""
function DroopControl(; params...)
    @named Pfil = Components.LowPassFilter(; τ=:τ_P, input=:P_meas, output=:P_fil)
    @named Pdroop = Components.Droop(; x_ref=:P_ref, K=:K_P, u_ref=:ω_ref, u=:ω, x=:P_fil)
    Psys = @connect Pfil.P_fil => Pdroop.P_fil outputs=:remaining name=:Psys

    @named Qfil = Components.LowPassFilter(; τ=:τ_Q, input=:Q_meas, output=:Q_fil)
    @named Qdroop = Components.Droop(; x_ref=:Q_ref, K=:K_Q, u_ref=:V_ref, u=:V, x=:Q_fil)
    Qsys = @connect Qfil.Q_fil => Qdroop.Q_fil outputs=:remaining name=:Qsys

    @named PQmeas = Components.Power(; u_r=:u_r_meas, u_i=:u_i_meas,
                                     i_r=:i_r_meas, i_i=:i_i_meas,
                                     P=:P_meas, Q=:Q_meas)

    refgen = Components.VRefGen()

    @named droopctrl = IOSystem(:autocon, [PQmeas, Psys, Qsys, refgen], outputs=:remaining)
    con = connect_system(droopctrl)

    return replace_vars(con, params)
end

function Synchronverter(; params...)
    # Frequency Loop
    @variables t ΔT(t) ω(t) θ(t)
    @parameters Te(t) P_ref Dp ω0 J
    dt = Differential(t)
    @named floop = IOBlock([ΔT ~ P_ref/ω0 - Dp*ω - Te,
                            dt(ω) ~ 1/J * ΔT,
                            dt(θ) ~ ω],
                           [Te],
                           [ω, θ])

    # Voltage Loop
    @variables t ΔQ(t) MfIf(t) V(t)
    @parameters u_r_meas(t) u_i_meas(t) Q(t) Q_ref Dq V_ref Kv
    @named vloop = IOBlock([V ~ sqrt(u_r_meas^2 + u_i_meas^2),
                            ΔQ ~ Q_ref - Q + Dq*(V_ref - V),
                            dt(MfIf) ~ ΔQ/Kv],
                           [Q, u_r_meas, u_i_meas],
                           [MfIf])
    BlockSystems.WARN[] = true

    # main model
    @variables t Te(t) u_r_ref(t) u_i_ref(t) Q(t)
    @parameters i_r_meas(t) i_i_meas(t) MfIf(t) ω(t) θ(t) ω0
    @named machine = IOBlock([Te      ~  sqrt(3/2) * MfIf * (cos(θ)*i_r_meas + sin(θ)*i_i_meas),
                              u_r_ref ~  sqrt(3/2) * MfIf * (ω0 + ω) * cos(θ),
                              u_i_ref ~  sqrt(3/2) * MfIf * (ω0 + ω) * sin(θ),
                              Q       ~ -sqrt(3/2) * MfIf * (ω0 + ω) * (-sin(θ)*i_r_meas + cos(θ)*i_i_meas)],
                             [ω, θ, i_r_meas, i_i_meas, MfIf],
                             [u_r_ref, u_i_ref, Q, Te])

    @named syncvert = IOSystem(:autocon, [floop, vloop, machine], outputs=:remaining, globalp=[:ω0])
    con = connect_system(syncvert)
    return replace_vars(con, params)
end

####
#### Inner Control Loops
####
"""
    PT1Source(; params...)

Create Schiffer Voltage source which follows angle directly but
amplitude with lag.
"""
function PT1Source(; params...)
    @variables t A(t) u_r(t) u_i(t)
    @parameters τ u_r_ref(t) u_i_ref(t)
    dt = Differential(t)
    @named PT1source = IOBlock([dt(A) ~ 1/τ*(√(u_r_ref^2 + u_i_ref^2) - A),
                                u_r ~ A/√(u_r_ref^2 + u_i_ref^2) * u_r_ref,
                                u_i ~ A/√(u_r_ref^2 + u_i_ref^2) * u_i_ref],
                               [u_r_ref, u_i_ref],
                               [u_r, u_i])

    if !isempty(params)
        PT1source = replace_vars(PT1source, params)
    end

    ui_meas = Components.UIMeas()

    @connect PT1source.(u_r, u_i) => ui_meas.(u_r, u_i) name=:PT1Source
end

"""
    PerfectSource()

Perfect Voltage source which follows the reference directly.
"""
function PerfectSource(; params...)
    @variables t u_r(t) u_i(t)
    @parameters u_r_ref(t) u_i_ref(t)
    @named Vsource = IOBlock([u_r ~ u_r_ref,
                              u_i ~ u_i_ref],
                             [u_r_ref, u_i_ref],
                             [u_r, u_i])

    ui_meas = Components.UIMeas()

    @connect Vsource.(u_r, u_i) => ui_meas.(u_r, u_i) name=:PT1Source
end
