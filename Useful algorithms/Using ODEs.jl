# Very important basis of how to use ODEs in JULIA!!
using ODE
using DifferentialEquations
using Plots

# Single variable
# function f(u, p, t)
#     return 1.01*u
# end
# u0 = 1 /2 
# tspan = (0,1)
# prob = ODEProblem(f, u0, tspan)
# sol = solve(prob, ode45(),reltol=1e-8,abstol=1e-8)
# plot(sol)


# Multiple variables
function f(u,p,t)
    x, v = u
    x_prime = v
    v_prime = -x
    return [x_prime, v_prime]
end
starty = [0.0,1.0]
time = (0.0, 4*pi)
prob = ODEProblem(f, starty, time)
sol = solve(prob, ode45(), reltol=1e-8,abstol=1e-30)
t = sol.t
u = sol.u[end]
# new_u = transpose(reduce(hcat, u))
# plot(t,new_u[:,1])
# plot!(t, new_u[:,2])


# Single variable looks like solve_ivp
# function g(t, u)
#     u_prime = 1.01 * u
#     return u_prime
# end

# u0 = 1/2
# time = 0:0.01:1
# t, y = ode45(g, u0, time)
# plot(sol)


# Multiple variables looks like solve_ivp
# function g(t, u)
#     x, v = u
#     x_prime = v
#     v_prime = -x
#     return [x_prime, v_prime]
# end

# starty = [0.0, 0.1]
# time = 0:0.1:4*pi
# t, y = ode45(g, starty, time)
# x = map(y -> y[1], y)
# v = map(y -> y[2], y)
# plot(t, x)
# plot!(t,v)


# Goes into negatives:
# function nega(y, p, t)
#     eps1 = p[1]
#     eps2 = p[2]
#     epsilon = eps1/eps2
#     return ((1-t)*y - y^2)/epsilon
# end

# starty = 1.0
# time= (0, 2)
# p = [1e-5, 10]
# prob = ODEProblem(nega, starty, time, p)
# sol = solve(prob, Tsit5(), isoutofdomain=(y,p,t)->any(x->x<0,y))
# plot(sol)

