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
        PT1source = set_p(PT1source, params)
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
