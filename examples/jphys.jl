using BlockSystems
using EMTSim
using Plots
using OrdinaryDiffEq
using NetworkDynamics
import EMTSim: Components


####
#### Schiffer voltage source
####
VS = EMTSim.SchifferSource()
@named P2C = Components.Polar2Cart(; mag=:u_mag, arg=:u_arg)
@named C2P = Components.Cart2Polar()

sys = @connect P2C.(x,y)=>VS.(û_r,û_i) VS.(u_r,u_i)=>C2P.(x,y) outputs=:remaining

@variables t
sys_tmp = set_input(sys, sys.u_mag => 1.0 + 0.1*sin(t))
sys_closed = set_input(sys_tmp, sys.u_arg => sin(2*t))
sys_f = set_p(sys_closed, :τ=>0.5)


f = ODEFunction(sys_f)
f.syms
u0 = Float64[1, 0, 1]
prob = ODEProblem(f, u0, (0, 5))
sol = solve(prob, Rosenbrock23())

plot(sol, vars=[:mag])
plot!(t->1 + 0.1*sin(t); label="u_mag")

plot(sol, vars=[:arg])
plot!(t->sin(2*t), label="u_arg")


####
#### Droop Control
####

# droop control with PT1 source
inner = EMTSim.PT1Source(τ=0.1)
outer = EMTSim.DroopControl(P_ref=1, Q_ref=0,
                            V_ref=1, ω_ref=0,
                            τ_P=0.1, K_P=1,
                            τ_Q=0.1, K_Q=1)

closed = ClosedLoop(inner, outer)
equations(closed)
odevert = ODEVertex(closed)

# droop control with perfect source
inner = EMTSim.PerfectSource()
outer = EMTSim.DroopControl(P_ref=1, Q_ref=0,
                            V_ref=1, ω_ref=0,
                            τ_P=0.1, K_P=1,
                            τ_Q=0.1, K_Q=1)

closed = ClosedLoop(inner, outer)
equations(closed)
odevert = ODEVertex(closed)

#

RMSPiLine(L=1, R=1, C1=1, C2=1)


StaticEdge

@variables t i_r_src(t) i_i_src(t) i_r_dst(t) i_i_dst(t)
# complex variables
@variables iC1(t) iL(t) iC2(t)
@parameters C1 C2 L R

a ~ conj(c)
c

@variables a b
c = a+im*b
conj(c)
imag(c)
real(c)
