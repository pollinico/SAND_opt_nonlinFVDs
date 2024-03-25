#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This code was written by Nicolò Pollini,                                %
# Technion - Israel Institute of Technology                               %  
#                                                                         %
#                                                                         %
# Contact: nicolo@technion.ac.il                                          %
#                                                                         %
# Code repository: https://github.com/pollinico/SAND_opt_nonlinFVDs       %
#                                                                         %
# Disclaimer:                                                             %
# The author reserves all rights but does not guarantee that the code is  %
# free from errors. Furthermore, the author shall not be liable in any    %
# event caused by the use of the program.                                 %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Import packages:
using JuMP, Ipopt, MAT, LinearAlgebra, Plots, Interpolations
using LaTeXStrings
gr()
#pyplot()

include("External_functions.jl")
include("Warm_start_ipopt.jl")
include("readIPOPToutput.jl")

WARM_START = true

# Read structural and loading input:
resultsFolder = "./res_2dof_Fy/"
folder = "./2dof_ex/"
Mfile = matread(folder*"M.mat")
Kfile = matread(folder*"K.mat")
Hfile = matread(folder*"H.mat")
Tfile = matread(folder*"T.mat")
efile = matread(folder*"e.mat")
LA02file = matread(folder*"LA02.mat")

M = Mfile["M"]
K = Kfile["K"]
H = Hfile["H"]
T = Tfile["T"]
e = efile["e"]
LA02 = LA02file["LA02"]
t   = LA02[1,:]
acc = LA02[2,:]

# Refine acc
itp = interpolate((t,), acc, Gridded(Linear()))
Δt_new = 0.006
t = collect(t[1] : Δt_new : t[end])
acc = itp(t)

tOld = t
tMax = 20.0 # seconds
t = t[tOld.<=tMax]
acc = acc[tOld.<=tMax]

plot(t, acc, α=1, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legend = false, margin=4Plots.mm)
xlabel!("t [s]")
ylabel!("a [m/s\$^2\$]")
savefig(resultsFolder*"LA02.png")

# Optimization problem settings:
β  = 1/4. # Newmark's parameters
γ  = 1/2.
Δt = t[2] - t[1]
dmax  = 15.0 # mm, max drift
dmaxVec = [13.0, 11.0, 9.0]

h = 3000.0 # mm
l = 5000.0 # mm
θ = asind(l/(h^2+l^2)^0.5)
cdMax = 5.0 # Max damping coefficient in kNs/mm
Fy = [169.0, 107.0] # Max structural force accepted
aParam = 5.0/100 # secondary stiffness post yield
nParam = 5.0 # shaerpness of transition elastic-inelastic

# Inherent damping matrix:
xi = 5.0/100
matA = K \ M
autoval = eigvals(matA)
omega2 = autoval.^-1
omega2 = sort!(omega2, rev=false)
omega = omega2.^.5
a0 = xi*2*omega[1]*omega[2]/(omega[1]+omega[2])
a1 = xi*2/(omega[1]+omega[2])
Cs = a0*M + a1*K

# Transform to displacements, from drfits:
m  = T' * M  * T / 1e3
cs = T' * Cs * T / 1e3
k  = T' * K  * T / 1e3
# I extract the diag terms of k in drift coordinates
diagK = diag(K) / 1e3

# Loading
P = - repeat(1e3*m*e,1,length(t)).*repeat(acc',2)

# Setting the problem dimension parameters:
ndof_d, ndof_u = size(H)
Nstep = length(t)
Ndamper = 2

# Scaling coefficient for problem equations
cEqEq     = 1.0e1
cFs       = 1.0e1
cStateVarU = 1.0e0
cStateVarV = 1.0e1
cStateVarA = 1.0e2

maxU = 3*dmax/cStateVarU # mm
maxUd = dmax/cStateVarU # mm
maxV =  Inf # 300.0/cStateVarV # mm/s
maxVd = Inf # 200.0/cStateVarV # mm/s
maxA = Inf # 10.0e3/cStateVarA # mm/s^2  
maxFs = Inf # 2*maximum(Fy)/cFs # kN
maxdFs = Inf # 5.0e3/cFs # maximum(diagK)
maxFd = Inf # 200.0/cFd # kN 

# Set up of the optimization problem in JuMP
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
set_optimizer_attribute(model, "linear_solver", "mumps") 
set_optimizer_attribute(model, "expect_infeasible_problem", "no")
set_optimizer_attribute(model, "max_iter", 500)
set_optimizer_attribute(model, "tol", 1e-6)
set_optimizer_attribute(model, "constr_viol_tol", 1e-6)
set_optimizer_attribute(model, "acceptable_tol", 5e-5)
set_optimizer_attribute(model, "acceptable_constr_viol_tol", 5e-5)
set_optimizer_attribute(model, "output_file", resultsFolder*"ipopt_out_Fy_dMax_"*string(Int(dmax))*".txt")
objHist = Float64[]
gHist   = Float64[]
iterIPOPT = Float64[]
timeIPOPT = Float64[]
Cd_solutions = zeros(2,0)
Obj_solutions = zeros(1,0)
MOI.set(model, Ipopt.CallbackFunction(), my_callback)

# Register function definitions
register(model, :sign, 1, sign; autodiff = true)

# Variables:
@variable(model, 1e-3   <= xd[1:Ndamper]         <= 1.0)
@variable(model, -maxU  <= xu[1:ndof_u,1:Nstep]  <= maxU)
@variable(model, -maxV  <= xv[1:ndof_u,1:Nstep]  <= maxV)
@variable(model, -maxA  <= xa[1:ndof_u,1:Nstep]  <= maxA)
if !WARM_START
    @variable(model, -maxUd <= d[1:ndof_d,1:Nstep] <= maxUd)
else
    d = @variable(model, [1:ndof_d,1:Nstep], lower_bound = -maxUd, upper_bound = maxUd)
end
@variable(model, -maxVd   <= vd[1:ndof_d,1:Nstep]     <= maxVd)
@variable(model, -maxFs   <= Fs[1:ndof_d,1:Nstep]     <= maxFs)
@variable(model, -maxdFs  <= dFs[1:ndof_d,1:Nstep-1]  <= maxdFs)
@variable(model, -maxFd   <= Fd[1:Ndamper,1:Nstep]    <= maxFd)
@variable(model, 0.0      <= τ                        <= Inf)

# Objective
@objective(model, Min, τ )

# Constraints:

for j in 1:Ndamper
    @NLconstraint(model, [i=2:Nstep], abs(Fd[j,i]) <= τ )
end

# Initial conditions
@constraint(model, u0,  xu[:,1] .== 0.0)
@constraint(model, v0,  xv[:,1] .== 0.0)
@constraint(model, Fs0, Fs[:,1] .== 0.0)
@constraint(model, Fd0, Fd[:,1] .== 0.0)
@constraint(model, a0,  xa[:,1] .== m \ (P[:,1]/cStateVarA - cs*zeros(ndof_u) - T'*Fd[:,1] - k*zeros(ndof_u)))

# Newmark-$\beta$
@constraint(model, New_v[i=1:Nstep-1], xv[:,i+1] .==  ( xv[:,i] + (1-γ) * Δt * xa[:,i]*cStateVarA/cStateVarV + γ * Δt * xa[:,i+1]*cStateVarA/cStateVarV) )
@constraint(model, New_u[i=1:Nstep-1], xu[:,i+1] .==  ( xu[:,i] + Δt * xv[:,i]*cStateVarV/cStateVarU + (Δt^2)*(1/2-β) * xa[:,i]*cStateVarA/cStateVarU + (Δt^2)*β * xa[:,i+1]*cStateVarA/cStateVarU) )

# Drift relation to displacments
@constraint(model, drift[i=1:Nstep], d[:,i] .== H * xu[:,i] )

# Drift velocity
@constraint(model, vdrift[i=1:Nstep], vd[:,i] .== H * xv[:,i] )

# Damper force
for j in 1:Ndamper
    @NLconstraint(model, [i=2:Nstep], Fd[j,i] ==  xd[j] * vd[j,i] )
end

# Incremental Fs, structural restoring force 
@constraint(model, f_struct[i=2:Nstep], Fs[:,i] .== (Fs[:,i-1] + dFs[:,i-1] * Δt ) )

# Stiffness matrix
for j in 1:ndof_d
    @NLconstraint(model, [i=1:Nstep-1],  dFs[j,i] == diagK[j] * (aParam + (1-aParam) * (1.0 - abs((Fs[j,i]*cFs-aParam*diagK[j]*d[j,i]*cStateVarU)/((1-aParam)*Fy[j]))^nParam * (0.5*sign(cStateVarV*(Fs[j,i]*cFs-aParam*diagK[j]*d[j,i]*cStateVarU)*vd[j,i]) + 0.5))  ) * vd[j,i] * cStateVarV / cFs  )
end

# Equilibrium equations
@constraint(model, eq_eq[i=2:Nstep], (m * xa[:,i] .* cStateVarA + cs * xv[:,i] .* cStateVarV + T' * Fd[:,i] * (cdMax*cStateVarV*sind(θ)^2) + T' * Fs[:,i] .* cFs) ./ cEqEq .== P[:,i] ./ cEqEq)

# Solution
optimize!(model)

# Save current solution
Cd_solutions = hcat(Cd_solutions, cdMax * value.(model[:xd]))
Obj_solutions = hcat(Obj_solutions, (cdMax*cStateVarV*sind(θ)) * objective_value(model))

# Extract data from output file of IPOPT
iter, time = readIpoptOutput(resultsFolder, "ipopt_out_Fy_dMax_"*string(Int(dmax))*".txt")
push!(iterIPOPT, iter)
push!(timeIPOPT, time)

# Use previous solution to start new optimization
if WARM_START
    for new_dmax in dmaxVec
        global iterIPOPT, timeIPOPT
        global iter, time
        global Cd_solutions, Obj_solutions
        set_optimal_start_values!(model)
        global dmax = new_dmax # I test the possibility of chaning the bounds and rerunning the model with a warm start from previous optimization for ipopt
        global maxUd = dmax/cStateVarU # mm/s
        for var_d in d
            set_lower_bound(var_d, -maxUd)
            set_upper_bound(var_d, +maxUd)
        end
        set_optimizer_attribute(model, "output_file", resultsFolder*"ipopt_out_Fy_dMax_"*string(Int(dmax))*".txt")
        optimize!(model)
        # Save current solution
        Cd_solutions = hcat(Cd_solutions, cdMax * value.(model[:xd]))
        Obj_solutions = hcat(Obj_solutions, (cdMax*cStateVarV*sind(θ)) * objective_value(model))
        # Extract data from output file of IPOPT
        iter, time = readIpoptOutput(resultsFolder, "ipopt_out_Fy_dMax_"*string(Int(dmax))*".txt")
        push!(iterIPOPT, iter)
        push!(timeIPOPT, time)
    end
end

solution_summary(model)

f = (cdMax*cStateVarV*sind(θ)) * objective_value(model)
Cd = cdMax * value.(model[:xd]) # [kNs/m]
u = cStateVarU * value.(model[:xu])
v = cStateVarV * value.(model[:xv])
a = cStateVarA * value.(model[:xa])
if WARM_START
    drift = cStateVarU * value.(d)
else
    drift = cStateVarU * value.(model[:d])
end
vdrift = cStateVarV * value.(model[:vd])
fs = cFs * value.(model[:Fs])
dfs = cFs*value.(model[:dFs])
fd = (cdMax*cStateVarV*sind(θ)) * value.(model[:Fd])

matwrite(resultsFolder*"u.mat", Dict("u" => u); compress = true)
matwrite(resultsFolder*"drift.mat", Dict("drift" => drift); compress = true)
matwrite(resultsFolder*"cd_opt.mat", Dict("cd_opt" => Cd_solutions); compress = true)
matwrite(resultsFolder*"obj_opt.mat", Dict("obj_opt" => Obj_solutions); compress = true)

println("")
println("Number of time steps: ", length(t))
println("Time step: ", Δt, " s")
println("Final objective function: ", round(f, digits=3), "  ")
println("Final cd1: ", round(Cd[1], digits=3), " kNs/mm")
println("Final cd2: ", round(Cd[2], digits=3), " kNs/mm")
println("Final max drift: ", round(maximum(abs.(drift)), digits=3), " mm")
println("Final max drift velocity: ", round(maximum(abs.(vdrift)), digits=3), " mm/s")
println("Final max displacement: ", round(maximum(abs.(u)), digits=3), " mm")
println("Final max velocity: ", round(maximum(abs.(v)), digits=3), " mm/s")
println("Final max acceleration: ", round(maximum(abs.(a)), digits=3), " mm/s^2")
println("Final max structural force: ", round(maximum(abs.(fs)), digits=3), " kN")
println("Final max dfs: ", round(maximum(abs.(dfs)), digits=3), " kN/s")
println("Final max damper force: ", round(maximum(abs.(fd)), digits=3), " kN")
println("Number of iterations in IPOPT: ", Int.(iterIPOPT))
println("Optimization time in IPOPT: ", round.(timeIPOPT, digits=1), " s")
println("Total number of iterations in IPOPT: ", sum(Int.(iterIPOPT)))
println("Total optimization time in IPOPT: ", round.(sum(timeIPOPT), digits=1), " s")

io = open(resultsFolder*"opt_sol.txt", "w")
write(io,string("Optimized solution for hysteretic structure with linear dampers","\n"))
write(io,string("Number of time steps: ", length(t),"\n"))
write(io,string("Time step: ", Δt, " s","\n"))
write(io,string("Final objective function: ", round(f, digits=3), "  ","\n"))
write(io,string("Final cd1: ", round(Cd[1], digits=3), " kNs/mm","\n"))
write(io,string("Final cd2: ", round(Cd[2], digits=3), " kNs/mm","\n"))
write(io,string("Final max drift: ", round(maximum(abs.(drift)), digits=3), " mm","\n"))
write(io,string("Final max drift velocity: ", round(maximum(abs.(vdrift)), digits=3), " mm/s","\n"))
write(io,string("Final max displacement: ", round(maximum(abs.(u)), digits=3), " mm","\n"))
write(io,string("Final max velocity: ", round(maximum(abs.(v)), digits=3), " mm/s","\n"))
write(io,string("Final max acceleration: ", round(maximum(abs.(a)), digits=3), " mm/s^2","\n"))
write(io,string("Final max structural force: ", round(maximum(abs.(fs)), digits=3), " kN","\n"))
write(io,string("Final max dfs: ", round(maximum(abs.(dfs)), digits=3), " kN/s","\n"))
write(io,string("Final max damper force: ", round(maximum(abs.(fd)), digits=3), " kN","\n"))
write(io,string("Number of iterations in IPOPT: ", Int.(iterIPOPT),"\n"))
write(io,string("Optimization time in IPOPT: ", round.(timeIPOPT, digits=1), " s","\n"))
write(io,string("Total number of iterations in IPOPT: ", sum(Int.(iterIPOPT)),"\n"))
write(io,string("Total optimization time in IPOPT: ", round.(sum(timeIPOPT), digits=1), " s","\n"))
close(io)

plot(t,  drift[1,:], label = "d1", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
plot!(t, drift[2,:], label = "d2", α=0.6, w=2)
plot!([t[1], t[end]], +[dmax, dmax], label="", linestyle = :dash, w=2, linecolor = :black)
plot!([t[1], t[end]], -[dmax, dmax], label="", linestyle = :dash, w=2, linecolor = :black)
xlabel!("t [s]")
ylabel!("d [mm]")
savefig(resultsFolder*"opt_drift_2D.png")

plot(t,  u[1,:], label = "u1", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
plot!(t, u[2,:], label = "u2", α=0.6, w=2)
xlabel!("t [s]")
ylabel!("u [mm]")
savefig(resultsFolder*"opt_u_2D.png")

plot(t,  v[1,:]/1.0e3, label = "v1", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
plot!(t, v[2,:]/1.0e3, label = "v2", α=0.6, w=2)
xlabel!("t [s]")
ylabel!("v [m/s]")
savefig(resultsFolder*"opt_v_2D.png")

plot(t,  a[1,:]/1.0e3, label = "a1", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
plot!(t, a[2,:]/1.0e3, label = "a2", α=0.6, w=2)
xlabel!("t [s]")
ylabel!("a [m/s\$^2\$]")
savefig(resultsFolder*"opt_a_2D.png")

plot(drift[1,:],  fs[1,:], label = "fs1", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
xlabel!("d [mm]")
ylabel!("fs1 [kN]")
savefig(resultsFolder*"opt_fs1_2D.png")

plot(drift[2,:],  fs[2,:], label = "fs2", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
xlabel!("d [mm]")
ylabel!("fs2 [kN]")
savefig(resultsFolder*"opt_fs2_2D.png")

plot(t,  fs[1,:], label = "fs1", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
plot!([t[1], t[end]], +[Fy[1], Fy[1]], label="", linestyle = :dash, w=2, linecolor = :black)
plot!([t[1], t[end]], -[Fy[1], Fy[1]], label="", linestyle = :dash, w=2, linecolor = :black)
xlabel!("t [s]")
ylabel!("fs1 [kN]")
savefig(resultsFolder*"opt_fs1_time_2D.png")

plot(t,  fs[2,:], label = "fs2", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
plot!([t[1], t[end]], +[Fy[2], Fy[2]], label="", linestyle = :dash, w=2, linecolor = :black)
plot!([t[1], t[end]], -[Fy[2], Fy[2]], label="", linestyle = :dash, w=2, linecolor = :black)
xlabel!("t [s]")
ylabel!("fs2 [kN]")
savefig(resultsFolder*"opt_fs2_time_2D.png")

plot(vdrift[1,:],  fd[1,:], label = "fd1", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
xlabel!("vd [mm/s]")
ylabel!("fd1 [kN]")
savefig(resultsFolder*"opt_fd1_vd_2D.png")

plot(vdrift[2,:],  fd[2,:], label = "fd2", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
xlabel!("vd [mm/s]")
ylabel!("fd2 [kN]")
savefig(resultsFolder*"opt_fd2_vd_2D.png")

plot(drift[1,:],  fd[1,:], label = "fd1", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
xlabel!("d [mm]")
ylabel!("fd1 [kN]")
savefig(resultsFolder*"opt_fd1_d_2D.png")

plot(drift[2,:],  fd[2,:], label = "fd2", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
xlabel!("d [mm]")
ylabel!("fd2 [kN]")
savefig(resultsFolder*"opt_fd2_d_2D.png")

plot(t, fd[1,:], label = "fd1", α=0.8, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
plot!(t, fd[2,:], label = "fd2", α=0.6, w=2)
plot!([minimum(t), maximum(t)], [f, f], label = "", w=2, linestyle=:dash, color=:black)
plot!([minimum(t), maximum(t)], [-f, -f], label = "", w=2, linestyle=:dash,  color=:black)
xlabel!("t [s]")
ylabel!("fd1, fd2 [kN]")
savefig(resultsFolder*"opt_fd12_t_2D.png")

plot(objHist, label = "", α=0.9, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
xlabel!("iter")
ylabel!("f")
savefig(resultsFolder*"opt_f_2D.png")

plot(gHist./maximum(abs.(gHist)), label = "", α=0.9, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm)
xlabel!("iter")
ylabel!("g")
savefig(resultsFolder*"opt_g_2D.png")

plot(Obj_solutions', label = "", α=0.9, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm, linestyle=:dash, markershape=:circle, markersize=6)
xlabel!("sub-problem")
ylabel!("τ [kN]")
savefig(resultsFolder*"opt_tau_2D.png")

plot(Cd_solutions[1,:], label = "cd1", α=0.9, w=2, size = (800, 400),  
xtickfont=15, ytickfont=15,  
guidefontsize=15, legendfont=15, margin=4Plots.mm, linestyle=:dash, markershape=:circle, markersize=6)
plot!(Cd_solutions[2,:], label = "cd2", α=0.6, w=2, linestyle=:dash, markershape=:circle, markersize=6)
xlabel!("sub-problem")
ylabel!("cd [kNs/mm]")
savefig(resultsFolder*"opt_cd_2D.png")
