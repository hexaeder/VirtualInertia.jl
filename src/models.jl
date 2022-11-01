export Slack, RMSPiLine, EMTRLLine, PT1Load, ConstPLoad, PT1PLoadEMT, SecondaryControlCS_PI, SecondaryControlCS_PT1, BSPiLine, ConstLoad

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
        StaticEdge(f=edgef, dim=4, coupling=:fiducial, sym=[:i_dst_r,:i_dst_i, :i_src_r, :i_src_i])
    end
end

function BSPiLine(; p_order=[], params...)
    @variables t i_src_r(t) i_src_i(t) i_dst_r(t) i_dst_i(t)
    @parameters u_src_r(t) u_src_i(t) u_dst_r(t) u_dst_i(t) R X B_src B_dst
    lineblock = IOBlock([i_src_r ~ real(-(u_src_r + im*u_src_i - u_dst_r - im*u_dst_i)/(R + im*X) - (im*B_src)*(u_src_r + im*u_src_i)),
                         i_src_i ~ imag(-(u_src_r + im*u_src_i - u_dst_r - im*u_dst_i)/(R + im*X) - (im*B_src)*(u_src_r + im*u_src_i)),
                         i_dst_r ~ real( (u_src_r + im*u_src_i - u_dst_r - im*u_dst_i)/(R + im*X) - (im*B_dst)*(u_dst_r + im*u_dst_i)),
                         i_dst_i ~ imag( (u_src_r + im*u_src_i - u_dst_r - im*u_dst_i)/(R + im*X) - (im*B_dst)*(u_dst_r + im*u_dst_i))],
                        [u_src_r, u_src_i, u_dst_r, u_dst_i],
                        [i_src_r, i_src_i, i_dst_r, i_dst_i];
                        name=:PiLine)
    lineblock = replace_vars(lineblock, params)
end

function EMTRLLine(; params...)
    @variables t i_r(t) i_i(t)
    @parameters u_src_r(t) u_src_i(t) u_dst_r(t) u_dst_i(t) R L ω0
    dt = Differential(t)
    lineblock = IOBlock([dt(i_r) ~  ω0 * i_i  - R/L * i_r + 1/L*(u_src_r - u_dst_r),
                         dt(i_i) ~ -ω0 * i_r  - R/L * i_i + 1/L*(u_src_i - u_dst_i)],
                        [u_src_r, u_src_i, u_dst_r, u_dst_i],
                        [i_r, i_i],
                        name=:RLLine)
    lineblock = replace_vars(lineblock, params)
    lineblock
end

function Slack()
    @variables t u_r(t) u_i(t)
    dt = Differential(t)
    slackblock = IOBlock([dt(u_r) ~ 0, dt(u_i) ~ 0], [], [u_r, u_i]; name=:slack)
end

function SwingEquation(; params...)
    @variables t ω(t) δ(t) P_el(t) u_r(t) u_i(t)
    @parameters i_r(t) i_i(t) P_ref D H ω0 u_mag
    dt = Differential(t)

    @named swing = IOBlock([P_el ~ i_i*u_i + i_r*u_r,
                            dt(ω) ~ ω0/(2*H) * (P_ref - P_el - D*ω),
                            dt(δ) ~ ω,
                            u_r ~ u_mag * cos(δ),
                            u_i ~ u_mag * sin(δ)],
                           [i_r, i_i], [u_r, u_i])
    swing = substitute_algebraic_states(swing)
    swing = replace_vars(swing, params)
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
    loadblock = replace_vars(loadblock, params)
end

function ConstLoad(;EMT=false, params...)
    cs = Components.PerfectCurrentSource(;name=:load_cs, P=:P_load, Q=:Q_load)
    bus = BusBar(cs; name=:pt1load, autopromote=true, EMT)
    bus = replace_vars(bus, params)
end

function PT1Load(;EMT=false, params...)
    cs = Components.PT1CurrentSource(;name=:load_cs, P=:P_load, Q=:Q_load)
    bus = BusBar(cs; name=:pt1load, autopromote=true, EMT)
    bus = replace_vars(bus, params)
end

function SecondaryControlCS_PI(; EMT=false, params...)
    cs = Components.PT1CurrentSource(P=:P_ref, Q=0)
    cs = make_input(cs, :P_ref)
    cs = replace_vars(cs; τ=0.001)

    pll = Components.ReducedPLL()
    pll = replace_vars(pll; ω_lp=1.32 * 2π*50, Kp=20.0, Ki=2.0)

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
    bus = replace_vars(bus, params)
end

function SecondaryControlCS_PT1(; EMT=false, params...)
    cs = Components.PT1CurrentSource(P=:P_ref_int, τ=:τ_cs, Q=0)
    cs = make_input(cs, :P_ref_int)

    lpf = Components.LowPassFilter(;:output=>:P_ref_int, :input=>:P_ref)
    lpf = make_iparam(lpf, :P_ref)

    cs = @connect lpf.P_ref_int => cs.P_ref_int outputs=:remaining name=:cs

    bus = BusBar(cs; autopromote=true, EMT)
    bus = replace_vars(bus, params)
end
