using EMTSim
using BlockSystems
using ModelingToolkit
using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using Plots
using Unitful
using LinearAlgebra
G
function get_rfcA(sym)
    @variables t x_d(t) x_q(t)
    @parameters x_r(t) x_i(t) δ(t)
    blk = IOBlock([x_d ~   x_r * cos(δ) + x_i * sin(δ),
                   x_q ~ - x_r * sin(δ) + x_i * cos(δ)],
                  [x_r, x_i, δ],
                  [x_d, x_q],
                  name=Symbol("rfcA_"*String(sym)))
    rename_vars(blk;
                x_r = Symbol(String(sym)*"_r"),
                x_i = Symbol(String(sym)*"_i"),
                x_d = Symbol(String(sym)*"_d"),
                x_q = Symbol(String(sym)*"_q"))
end

function get_rfcB(sym)
    @variables t x_r(t) x_i(t)
    @parameters x_d(t) x_q(t) δ(t)
    blk = IOBlock([x_r ~ x_d * cos(δ) - x_q * sin(δ),
                   x_i ~ x_d * sin(δ) + x_q * cos(δ)],
                  [x_d, x_q, δ],
                  [x_r, x_i],
                  name=Symbol("rfcB_"*String(sym)))
    rename_vars(blk;
                x_r = Symbol(String(sym)*"_r"),
                x_i = Symbol(String(sym)*"_i"),
                x_d = Symbol(String(sym)*"_d"),
                x_q = Symbol(String(sym)*"_q"))
end

####
#### LCL Filter
####
@variables t i_f_r(t) i_f_i(t) i_g_r(t) i_g_i(t) V_C_r(t) V_C_i(t)
@parameters ω0 Rf Rg Lf Lg C V_I_r(t) V_I_i(t) V_g_r(t) V_g_i(t)
dt = Differential(t)

LCL = IOBlock([dt(i_f_r) ~ -Rf/Lf*i_f_r + ω0*i_f_i + 1/Lf*(-V_C_r + V_I_r),
               dt(i_f_i) ~ -Rf/Lf*i_f_i - ω0*i_f_r + 1/Lf*(-V_C_i + V_I_i),
               dt(i_g_r) ~ -Rg/Lg*i_g_r + ω0*i_g_i + 1/Lg*(V_C_r - V_g_r),
               dt(i_g_i) ~ -Rg/Lg*i_g_i - ω0*i_g_r + 1/Lg*(V_C_i - V_g_i),
               dt(V_C_r) ~ 1/C * (i_f_r - i_g_r) + ω0*V_C_i,
               dt(V_C_i) ~ 1/C * (i_f_i - i_g_i) - ω0*V_C_r],
              [V_g_r, V_g_i, V_I_r, V_I_i],
              [i_f_r, i_f_i, i_g_r, i_g_i, V_C_r, V_C_i],
              name = :LCL)

####
#### Current Controller
####
@variables t γ_d(t) γ_q(t) V_I_d(t) V_I_q(t)
@parameters KP KI i_f_d(t) i_f_q(t) i_f_ref_d(t) i_f_ref_q(t) Lf
dt = Differential(t)

CC = IOBlock([dt(γ_d) ~ i_f_ref_d - i_f_d,
              dt(γ_q) ~ i_f_ref_q - i_f_q,
              V_I_d ~  Lf*ω0*i_f_q + KP*(i_f_ref_d - i_f_d) + KI*γ_d,
              V_I_q ~ -Lf*ω0*i_f_d + KP*(i_f_ref_q - i_f_q) + KI*γ_q],
             [i_f_d, i_f_q, i_f_ref_d, i_f_ref_q],
             [V_I_d, V_I_q],
             name=:CC)

@variables t γ_d(t) γ_q(t) i_f_ref_d(t)  i_f_ref_q(t)
@parameters KP KI F V_C_d(t) V_C_q(t) V_C_ref_d(t) V_C_ref_q(t) i_g_d(t) i_g_q(t) ω0 C
dt = Differential(t)

####
#### Voltage Controller
####
VC = IOBlock([dt(γ_d) ~ V_C_ref_d - V_C_d,
              dt(γ_q) ~ V_C_ref_q - V_C_q,
              i_f_ref_d ~ F*i_g_d + C*ω0*V_C_q + KP*(V_C_ref_d - V_C_d) + KI*γ_d,
              i_f_ref_q ~ F*i_g_q - C*ω0*V_C_d + KP*(V_C_ref_q - V_C_q) + KI*γ_q],
             [i_g_d, i_g_q, V_C_d, V_C_q, V_C_ref_d, V_C_ref_q],
             [i_f_ref_d, i_f_ref_q],
             name=:VC)


####
#### RFC Blocks
####
rfc_i_f = get_rfcA(:i_f)
rfc_i_g = get_rfcA(:i_g)
rfc_V_C = get_rfcA(:V_C)
rfc_V_I = get_rfcB(:V_I)

####
#### System Creation
####
system = IOSystem(:autocon,
                  [LCL,CC, VC, rfc_i_f, rfc_i_g, rfc_V_C, rfc_V_I],
                  outputs=[LCL.i_g_r, LCL.i_g_i],
                  globalp=[:δ, :C, :ω0, :Rf, :Rg, :Lf, :Lg],
                  name=:inv)

sys = connect_system(system)
sys = set_input(sys, :δ=>0)
# sys = set_input(sys, :V_C_ref_d=>1)
# sys = set_input(sys, :V_C_ref_q=>0)

params = Dict(
    :VC₊KI => 1,
    :CC₊KI => 1,
    :VC₊KP => 1,
    :CC₊KP => 1,
    :F => 1,
    :Rf => 0.01,
    :Rg => 0.01,
    :Lf => 10e-6,
    :Lg => 10e-6,
    :C => 1e-6,
    :ω0 => 2π*50,
)

sysp = set_p(sys, params) |> simplify

A, B, C, D = identify_lti(sysp)
A = fixdiv(A)


ss = StateSpace(A, B, Matrix(I,10,10), zeros(10,4))

OG = gram(ss,:o)
isposdef(OG)
CG = gram(ss,:c)
isposdef(CG)

sysbt, G, T = baltrunc(ss)
bodeplot(ss)

pzmap(ss)

eigen(A)


sysclosed = set_input(sysp, :V_g_r=>0.8)
sysclosed = set_input(sysclosed, :V_g_i=>0.0)

f = ODEFunction(sysclosed)
u0 = zeros(length(f.syms))
prob = ODEProblem(f, u0, (0, 10.))
sol = solve(prob, Rodas5())

f.syms

plot(sol, vars=[:i_g_r, :i_g_i])
# plot!(sol, vars=[:i_f_r, :i_f_i])
plot(sol, vars=[:V_C_r, :V_C_i])
plot(sol, vars=[:V_C_r, :V_C_i])
plot(sol, vars=[:CC₊γ_d, :CC₊γ_q])
plot(sol, vars=[:VC₊γ_d, :VC₊γ_q])

xlims!(0, 0.001)

f.syms
eigen(A)

Ctrb = ctrb(A, B)
rank(Ctrb)

det(Ctrb)

Obsv = obsv(A, I)
rank(Obsv)
Obsv

rank(Obsv[1:10,:])
svd(Obsv[1:10,:])

svd(Ctrb)

rank(Obsv)

rank(Obsv)

Matrix(I, 4, 4)
