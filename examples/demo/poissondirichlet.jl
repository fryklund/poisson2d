using bieps2d
using FMMLIB2D
using LinearAlgebra
import LinearMaps
import IterativeSolvers
#using Plots
using SpecialFunctions
include("../../src/julia/laplace.jl")
include("../../src/julia/helsingquad_laplace.jl")

eulergamma = Base.MathConstants.eulergamma


Lgrid = 2.0
sigma = 0.25
Ngrid = 642

# Poisson problem
# Exact solution from Gaussian blob rhs (from Vico & Greengard)

#     x0 = 0.123
#    y0 = 0.5193
#   r(x,y) = sqrt( (x-x0)^2 + (y-y0)^2 )
#  ufuncr(r) = -1/(4.0*pi)*(expint(r^2/(2.0*sigma^2)) + log(r^2))
# ulim0 = (eulergamma + log(1/(2*sigma^2)))/(4*pi)
# ffuncr(r) = -exp(-r^2/(2*sigma^2))/(2*pi*sigma^2)    
# ufunc(x, y) = ufuncr(r(x,y))
# ffunc(x, y) = ffuncr(r(x,y))

#= Volume FMM example
rsig = 1.0/1000.0

#     Gaussian source centers and variance
dpars = Vector{Float64}(undef,15)
dpars[1] = 0.1
dpars[2] = 0.15

dpars[3] = rsig

dpars[4] = -0.18
dpars[5] = 0.0

dpars[6] = rsig/2.1

dpars[7] = 0.178
dpars[8] = -0.1

dpars[9] = rsig/4.5

dpars[10] = -0.112
dpars[11] = 0.2

dpars[12] = rsig/1.2
ufunc(x,y) = exp(-((x-dpars[1])^2 + (y-dpars[2])^2)/dpars[3]) + exp(-((x-dpars[4])^2 + (y-dpars[5])^2)/dpars[6]) + exp(-((x-dpars[7])^2 + (y-dpars[8])^2)/dpars[9]) + exp(-((x-dpars[10])^2 + (y-dpars[11])^2)/dpars[12])

ffunc(x,y) = exp(-((x-dpars[1])^2 + (y-dpars[2])^2)/dpars[3])*(((x-dpars[1])^2 + (y-dpars[2])^2)/dpars[3] - 1)*4/dpars[3] + exp(-((x-dpars[4])^2 + (y-dpars[5])^2)/dpars[6])*(((x-dpars[4])^2 + (y-dpars[5])^2)/dpars[6] - 1)*4/dpars[6] + exp(-((x-dpars[7])^2 + (y-dpars[8])^2)/dpars[9])*(((x-dpars[7])^2 + (y-dpars[8])^2)/dpars[9] - 1)*4/dpars[9] + exp(-((x-dpars[10])^2 + (y-dpars[11])^2)/dpars[12])*(((x-dpars[10])^2 + (y-dpars[11])^2)/dpars[12] - 1)*4/dpars[12]

=#

L = Lgrid * 2

#ufunc(x,y) = 1/(cos((2*pi*(x + y))/L) + 2)
#ffunc(x,y) = -(4*pi^2*(-3 - 4*cos((2*pi*(x + y))/L) + cos((4*pi*(x + y))/L)))/(L^2*(2 + cos((2*pi*(x + y))/L))^3)

k = 40*pi #Bruno
ufunc(x,y) = sin(k*x)*sin(k*y) / (2 * k^2)
ffunc(x,y) = -sin(k*x)*sin(k*y) #* (2 * k^2)

# Discretize
numpanels = 200
panelorder = 16
#    curve = AnalyticDomains.starfish(n_arms = 3, amplitude=0.3)
curve = AnalyticDomains.kite()
dcurve = CurveDiscretization.discretize(curve, numpanels, panelorder)
arclength = sum(dcurve.dS)

# Volume grid
xgrid = range(-Lgrid, Lgrid, length=Ngrid+1)
ygrid = range(-Lgrid, Lgrid, length=Ngrid+1)
xgrid = xgrid[1:end-1]
ygrid = ygrid[1:end-1]
X, Y = PUX.ndgrid(xgrid, ygrid)
xe = [vec(X) vec(Y)]
xet = copy(xe')
interior, interior_near = CurveDiscretization.interior_points(dcurve, xet)
# Setup input and reference in grid
FIN = ffunc.(X, Y)
UREF = ufunc.(X, Y)

Fedges=[FIN[1,:];FIN[end,:];FIN[:,1];FIN[:,end]]
Fgridresolution = norm(Fedges, Inf) / norm(vec(FIN), Inf)
@show Fgridresolution

# PUX
pux_ep = 2
pux_P = 36
# Prepare
PUXStor = PUX.pux_params(Lgrid, Ngrid, pux_ep, pux_P, arclength)
# Get equispaced points
_, z_equi = CurveDiscretization.traparclen(curve, PUXStor.numcenters)
z_grid = dcurve.points[1,:] + 1im*dcurve.points[2,:]    
# Precompute
PUX.pux_setup!(xe, interior, z_grid, z_equi, PUXStor)


# Extend
println(" * PUX")
@time FEXT = PUX.pux_eval(FIN, [PUXStor])

# Poisson solve
println(" * Poisson solve")
xt, yt = dcurve.points[1,:], dcurve.points[2,:]
@time Up, Up_bdry = FreeSpace.fs_poisson(FEXT, 2*Lgrid, xt, yt)

# Laplace solve
bdry_cond = ufunc.(xt, yt)
mod_bdry_cond = bdry_cond - Up_bdry
@show norm(mod_bdry_cond)

println(" * Integral equation solve")
S = zeros(Float64,2,0)
flhs = system_matvec(dcurve,S)
LHS = LinearMaps.LinearMap(flhs, dcurve.numpoints)
rhs = mod_bdry_cond
sol, gmlog = IterativeSolvers.gmres(LHS, rhs; reltol=eps(), log=true)
gmresidual = norm(LHS*sol-rhs, Inf)
@show gmresidual



density = sol

println(" * Layer pot eval")
@time uh = layer_potential_fmm(dcurve, density, xet[:, interior],interior_near,S; slp = 0, dlp = 1)

Uh = zeros(size(X))
Uh[interior] = uh
U = Uh + Up
E = U - UREF
E[.!interior] .= 0
U[.!interior] .= 0
UREF[.!interior] .= 0    
Einf = norm(E[interior], Inf)
h = X[2,1]-X[1,1]
EL2 = sqrt(sum(E[interior].^2)*h^2)
@show Einf
@show EL2

