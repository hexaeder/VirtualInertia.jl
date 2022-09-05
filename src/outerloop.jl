"""
    DroopControl(; params...)

Return block for droop control outer ocntrol.
"""
function DroopControl(; params...)
    @named Pfil = Components.LowPassFilter(; τ=:τ_P, input=:P, output=:P_fil)
    @named Pdroop = Components.Droop(; x_ref=:P_ref, K=:K_P, u_ref=:ω_ref, u=:ω, x=:P_fil)
    Psys = @connect Pfil.P_fil => Pdroop.P_fil outputs=:remaining name=:Psys

    @named Qfil = Components.LowPassFilter(; τ=:τ_Q, input=:Q, output=:Q_fil)
    @named Qdroop = Components.Droop(; x_ref=:Q_ref, K=:K_Q, u_ref=:V_ref, u=:V, x=:Q_fil)
    Qsys = @connect Qfil.Q_fil => Qdroop.Q_fil outputs=:remaining name=:Qsys

    @named PQmeas = Components.Power(; u_r=:u_r_meas, u_i=:u_i_meas,
                                     i_r=:i_r_meas, i_i=:i_i_meas)

    refgen = Components.VRefGen()

    @named droopctrl = IOSystem(:autocon, [PQmeas, Psys, Qsys, refgen], outputs=:remaining)
    con = connect_system(droopctrl)

    if !isempty(params)
        con = set_p(con, params)
    end

    return con
end
