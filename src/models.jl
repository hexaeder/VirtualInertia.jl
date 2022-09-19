export Slack, RMSPiLine, EMTRLLine, PT1PLoad, ConstPLoad, PT1PLoadEMT, SecondaryControlCS

function RMSPiLine(; L, R, C1, C2)
    let ω=2π*50, L=L, R=R, C1=C1, C2=C2
        edgef = function(u, src, dst, p, t)
            Vsrc = Complex(src[1], src[2])
            Vdst = Complex(dst[1], dst[2])
            imain = (Vsrc - Vdst)/(R + im*ω*L)
            iC1 = Vsrc*im*ω*C1
            iC2 = Vdst*im*ω*C2
            isrc = - (imain + iC1) # positiv current means flow to node
            idst = imain - iC2
            u[1] = real(idst)
            u[2] = imag(idst)
            u[3] = real(isrc)
            u[4] = imag(isrc)
        end
        StaticEdge(f=edgef, dim=4, coupling=:fiducial, sym=[:i_r_dst,:i_i_dst, :i_r_src, :i_i_src])
    end
end

function EMTRLLine(; params...)
    @variables t i_r(t) i_i(t)
    @parameters u_r_src(t) u_i_src(t) u_r_dst(t) u_i_dst(t) R L ω0
    dt = Differential(t)
    lineblock = IOBlock([dt(i_r) ~  ω0 * i_i  - R/L * i_r + 1/L*(u_r_src - u_r_dst),
                         dt(i_i) ~ -ω0 * i_r  - R/L * i_i + 1/L*(u_i_src - u_i_dst)],
                        [u_r_src, u_i_src, u_r_dst, u_i_dst],
                        [i_r, i_i],
                        name=:RLLine)
    lineblock = set_p(lineblock, params)
    if !isempty(lineblock.iparams)
        @warn "There are open parameters on this line: $(lineblock.iparams)"
    end
    ODEEdge(lineblock)
end

function Slack()
    @variables t u_r(t) u_i(t)
    dt = Differential(t)
    slackblock = IOBlock([dt(u_r) ~ 0, dt(u_i) ~ 0], [], [u_r, u_i]; name=:slack)
    ODEVertex(slackblock)
end

function PT1PLoad(;params...)
    @variables t i_inj_r(t) i_inj_i(t) P_meas(t) Q_meas(t) u_r(t) u_i(t)
    @parameters P_ref τ i_r(t) i_i(t)
    dt = Differential(t)
    loadblock = IOBlock([P_meas ~ u_r*i_inj_r + u_i*i_inj_i,
                         Q_meas ~ i_inj_r*u_i - i_inj_i*u_r,
                         dt(i_inj_r) ~ 1/τ * (P_ref * u_r/(u_r^2 + u_i^2) - i_inj_r),
                         dt(i_inj_i) ~ 1/τ * (P_ref * u_i/(u_r^2 + u_i^2) - i_inj_i),
                         0 ~ i_inj_r + i_r,
                         0 ~ i_inj_i + i_i],
                        [i_r, i_i], [u_r, u_i],
                        name=:load)
    loadblock = substitute_algebraic_states(loadblock)
    loadblock = set_p(loadblock, params)

    ODEVertex(loadblock, ModelingToolkit.getname.(loadblock.iparams))
end

function ConstPLoad(;params...)
    @variables t P_meas(t) Q_meas(t) u_r(t) u_i(t)
    @parameters P_ref i_r(t) i_i(t)
    loadblock = IOBlock([P_meas ~ -u_r*i_r - u_i*i_i,
                         Q_meas ~ -i_r*u_i + i_i*u_r,
                         0 ~ Q_meas,
                         0 ~ P_ref - P_meas],
                        [i_r, i_i], [u_r, u_i],
                        name=:load)
    loadblock = substitute_algebraic_states(loadblock)
    loadblock = set_p(loadblock, params)

    ODEVertex(loadblock, ModelingToolkit.getname.(loadblock.iparams))
end

function PT1PLoadEMT(;params...)
    @variables t i_r(t) i_i(t) P_meas(t) Q_meas(t)
    @parameters u_r(t) u_i(t) P_ref τ
    dt = Differential(t)
    loadblock = IOBlock([P_meas ~ u_r*i_r + u_i*i_i,
                         Q_meas ~ i_r*u_i - i_i*u_r,
                         dt(i_r) ~ 1/τ*(P_ref * u_r/(u_r^2 + u_i^2) - i_r),
                         dt(i_i) ~ 1/τ*(P_ref * u_i/(u_r^2 + u_i^2) - i_i)],
                        [u_r, u_i], [i_r, i_i],
                        name=:load)
    loadblock = substitute_algebraic_states(loadblock)
    bus = BusBar(loadblock; name=:loadbus, autopromote=true)
    bus = set_p(bus, params)
    ODEVertex(bus, ModelingToolkit.getname.(bus.iparams))
end

function SecondaryControlCS(; EMT=false, params...)
    cs = Components.PT1CurrentSource()
    cs = set_p(cs; τ=0.01)

    pll = Components.ReducedPLL()
    pll = set_p(pll; ω_lp=1.32 * 2π*50, Kp=20.0, Ki=2.0)

    @variables t P_ref_pi(t) μ(t)
    @parameters P_ref ω(t) Ki Kp
    dt = Differential(t)
    ctrl = IOBlock([dt(μ) ~ -ω,
                    P_ref_pi ~ P_ref + Ki*μ - Kp*ω],
                   [ω], [P_ref_pi],
                   name=:PrefICtrl)

    @named pisrc = IOSystem([pll.ω_pll=>ctrl.ω,
                             ctrl.P_ref_pi=>cs.P_ref],
                            [pll, cs, ctrl],
                            outputs=[cs.i_r, cs.i_i],
                            globalp=[:u_r, :u_i])
    pisrc = connect_system(pisrc)

    bus = BusBar(pisrc; autopromote=true, EMT)
    bus = set_p(bus, params)
    ODEVertex(bus, ModelingToolkit.getname.(bus.iparams))
end
