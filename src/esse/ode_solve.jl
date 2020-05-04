#=

Create ODE numerical integration function

Ignacio Quintero Mächler

t(-_-t)

September 19 2017

=#





"""
    solvef(int::OrdinaryDiffEq.ODEIntegrator,
           u  ::Array{Float64,1},
           ti ::Float64,
           tf ::Float64)

Solve an IDE integrator for new `u` and and `ti` and `tf`.
"""
function solvef(int::OrdinaryDiffEq.ODEIntegrator,
                u  ::Array{Float64,1},
                ti ::Float64,
                tf ::Float64)
  @inbounds begin
    reinit!(int, u, t0 = ti, tf = tf)
    return solve!(int).u[1]
  end
end














# odef = ode_fun


# u0 = copy(X[1])
# p0 = copy(p)

# prob = ODEProblem(odef0, u0, (0.0,1.0), p0)
# int = init(prob,
#            Tsit5(),
#            save_everystep  = false, 
#            calck           = false,
#            force_dtmin     = true,
#            save_start      = false,
#            initialize_save = false,
#            maxiters        = 100_000_000,
#            verbose         = false)
# solve!(int).u[1]

# @benchmark begin 
#   int.p = p
#   reinit!(int, u0, t0 = 0.0, tf = 10.0)
#   solve!(int)
# end


# @benchmark begin 
#   reinit!(int, u0, t0 = 0.0, tf = 10.0)
#   solve!(int)
# end




# uS = SVector{6}(u0)
# p0 = copy(p)

# probS = ODEProblem(odefSS, uS, (0.0,1.0), p0)
# intS = init(prob,
#            Tsit5(),
#            save_everystep  = false, 
#            calck           = false,
#            force_dtmin     = true,
#            save_start      = false,
#            initialize_save = false,
#            maxiters        = 100_000_000,
#            verbose         = false)
# solve!(intS).u[1]

# @benchmark begin 
#   #copyto!(intS.p,p)
#   reinit!(intS, uS, t0 = 0.0, tf = 10.0)
#   solve!(intS)
# end

# @benchmark intS.p = p

# @benchmark copyto!(intS.p, p)
# @benchmark unsafe_copyto!(intS.p, 1, p, 1, 11)
# @benchmark intS.p = p


# p[1] = 1.0






# function odef0(du, u, p, t)
#   @inbounds begin
#     #af!(t, r)
#     #eaft[1] = p[1] * exp(p[10] * 1.0)
#     #eaft[2] = p[2] * exp(p[11] * 2.0)
#     du[1] = -1.0 * (p[1] + p[4] + p[6]) * u[1] + +(p[6] * u[3]) + 0.0 + 2.0 * (p[1] * u[4] * u[1])
#     du[2] = -1.0 * (p[2] + p[5] + p[7]) * u[2] + +(p[7] * u[3]) + 0.0 + 2.0 * (p[2] * u[5] * u[2])
#     du[3] = -1.0 * (1.0 * p[3] + p[2] + p[9] + p[1] + p[8]) * u[3] + (p[9] * u[1] + p[8] * u[2]) + 0.0 + (p[2] * (u[5] * u[3] + u[6] * u[2]) + p[1] * (u[4] * u[3] + u[6] * u[1])) + p[3] * (u[4] * u[2] + u[5] * u[1])
#     du[4] = -1.0 * (p[1] + p[4] + p[6]) * u[4] + +(p[4]) + +(p[6] * u[6]) + 0.0 + p[1] * u[4] ^ 2
#     du[5] = -1.0 * (p[2] + p[5] + p[7]) * u[5] + +(p[5]) + +(p[7] * u[6]) + 0.0 + p[2] * u[5] ^ 2
#     du[6] = -1.0 * (1.0 * p[3] + p[2] + p[9] + p[1] + p[8]) * u[6] + (p[9] * u[4] + p[8] * u[5]) + 0.0 + (p[2] * u[5] * u[6] + p[1] * u[4] * u[6]) + p[3] * +(u[4] * u[5])
#     return nothing
#   end
# end


# function odefSS(u, p, t)
#   @inbounds begin
#     #af!(t, r)
#     #eaft[1] = p[1] * exp(p[10] * 1.0)
#     #eaft[2] = p[2] * exp(p[11] * 2.0)
#     du1 = -1.0 * (p[1] + p[4] + p[6]) * u[1] + +(p[6] * u[3]) + 0.0 + 2.0 * (p[1] * u[4] * u[1])
#     du2 = -1.0 * (p[2] + p[5] + p[7]) * u[2] + +(p[7] * u[3]) + 0.0 + 2.0 * (p[2] * u[5] * u[2])
#     du3 = -1.0 * (1.0 * p[3] + p[2] + p[9] + p[1] + p[8]) * u[3] + (p[9] * u[1] + p[8] * u[2]) + 0.0 + (p[2] * (u[5] * u[3] + u[6] * u[2]) + p[1] * (u[4] * u[3] + u[6] * u[1])) + p[3] * (u[4] * u[2] + u[5] * u[1])
#     du4 = -1.0 * (p[1] + p[4] + p[6]) * u[4] + +(p[4]) + +(p[6] * u[6]) + 0.0 + p[1] * u[4] ^ 2
#     du5 = -1.0 * (p[2] + p[5] + p[7]) * u[5] + +(p[5]) + +(p[7] * u[6]) + 0.0 + p[2] * u[5] ^ 2
#     du6 = -1.0 * (1.0 * p[3] + p[2] + p[9] + p[1] + p[8]) * u[6] + (p[9] * u[4] + p[8] * u[5]) + 0.0 + (p[2] * u[5] * u[6] + p[1] * u[4] * u[6]) + p[3] * +(u[4] * u[5])
#   end

#   @SVector [du1,du2,du3,du4,du5,du6]
# end




