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



Ns = 20
MXM = 40 #6^Ns
psi = makePsi0(Ns)

Id = Matrix{Float64}(I, 3, 3)
base = Array{Float64,2}[Id for i = 1:Ns]



#Tuned to AKLT point without constant term expect energy  of -(Ns-1)/3. Very efficient construction, pre-allocates memory.
@time H = XVBSmake(Ns,(1/2),(1/6))






nD = 25 # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,H,maxm=MXM,sweeps=500,cutoff=1E-9, storeD=dummy_D,allSvNbond=true)
popat!(params.SvNvec,1) #this is here because SvNnev always produces a duplicated first entry
