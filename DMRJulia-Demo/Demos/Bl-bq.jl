#Now let's move to a slightly more sophisticated example
#########################################################################
#
#  Density Matrix Renormalization Group (and other methods) in julia (DMRjulia)
#                              v0.1
#
#########################################################################
# Made by Thomas E. Baker (2018)
# See accompanying license with this program
# This code is native to the julia programming language (v1.1.0) or (v1.5)
#

#=
IF RUNNING ON CLUSTER UNCOMMENT!
import Pkg
Pkg.add("LinearAlgebra")
Pkg.add("DelimitedFiles")
=#
using LinearAlgebra
using DelimitedFiles



path = "../"
include(join([path,"DMRjulia.jl"]))
using .DMRjulia
using LinearAlgebra
using Plots
using LaTeXStrings


Ns = 20
MXM = 40 #6^Ns
psi = makePsi0(Ns)

Id = Matrix{Float64}(I, 3, 3)
base = Array{Float64,2}[Id for i = 1:Ns]



############################################################################################################
θ = atan(1/3)
#Tuned to AKLT point without constant term, rescaled from previous.
@time H = XVBSmake(Ns,cos(θ),sin(θ))

nD = 25 # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,H,maxm=MXM,sweeps=500,cutoff=1E-9, storeD=dummy_D,allSvNbond=true)
popat!(params.SvNvec,1) #this is here because SvNnev always produces a duplicated first entry

plot(2:Ns,params.SvNvec, label = false, ylabel = "von Neumann Entanglement Entropy", xlabel = "Bipartition placement (between site i and i+1)")



############################################################################################################
#now let's tune away from the AKLT point, starting with the ferromagnetic point. The model is purely negative bilinear coupling so spins should be aligned with one another and the magnitude of \vec{S}_{i}\cdot\vec{S}_{i+1} should be as large as possible. This is a so-called gapless phase-- one where the energy gap between the gs and first excited state decays to zero in the thermodynamic limit. Expect a dome shaped entanglement structure.
θ = Float64(pi)*(0.875)
#Tuned to AKLT point without constant term, rescaled from previous.
@time H = XVBSmake(Ns,cos(θ),sin(θ))
psi = makePsi0(Ns)
nD = 25 # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,H,maxm=MXM,sweeps=500,cutoff=1E-9, storeD=dummy_D,allSvNbond=true)
popat!(params.SvNvec,1) #this is here because SvNnev always produces a duplicated first entry

plot(2:Ns,params.SvNvec, label = false, ylabel = "von Neumann Entanglement Entropy", xlabel = "Bipartition placement (between site i and i+1)")

#let's tune to the dimerized phase. We expect a dimerized ordering (ie singlet formations between nearest neighbour sites). This will lead to a jagged entanglement structure as we cleave across dimerized sites and undimerized sites with out bipartition.
############################################################################################################
θ = -Float64(pi)/2
@time H = XVBSmake(Ns,cos(θ),sin(θ))
psi = makePsi0(Ns)
nD = 25 # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,H,maxm=MXM,sweeps=40,cutoff=1E-9, storeD=dummy_D,allSvNbond=true)
popat!(params.SvNvec,1) #this is here because SvNnev always produces a duplicated first entry

plot(2:Ns,params.SvNvec, label = false, ylabel = "von Neumann Entanglement Entropy", xlabel = "Bipartition placement (between site i and i+1)")



#let's tune to the trimerized phase. This is an example of gapless phase. Expect similar entanglement structure to ferromagnetic phase.
############################################################################################################
θ = Float64(pi)*(0.75)
@time H = XVBSmake(Ns,cos(θ),sin(θ))
psi = makePsi0(Ns)
MXM = 65
nD = 100 # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,H,maxm=MXM,sweeps=50,cutoff=1E-6, storeD=dummy_D,allSvNbond=true)
popat!(params.SvNvec,1) #this is here because SvNnev always produces a duplicated first entry

plot(2:Ns,params.SvNvec, label = false, ylabel = "von Neumann Entanglement Entropy", xlabel = "Bipartition placement (between site i and i+1)")
