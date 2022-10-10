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

    return set_p(con, params)
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
    return set_p(con, params)
end
