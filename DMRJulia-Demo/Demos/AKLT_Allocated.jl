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



#Tuned to AKLT point without constant term expect energy  of -(Ns-1)/3. Very efficient construction, pre-allocates memory.
@time H = XVBSmake(Ns,(1/2),(1/6))






nD = 25 # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,H,maxm=MXM,sweeps=500,cutoff=1E-9, storeD=dummy_D,allSvNbond=true)
popat!(params.SvNvec,1) #this is here because SvNnev always produces a duplicated first entry




#it is very common to use .txt files to save data that outputs from numerical studies.
savearray = [Vector(2:Ns),params.SvNvec]
filename = string("AKLT_EE_",Ns,"_SITE.txt")
touch(filename)
open(filename,"w") do io
    writedlm(filename, savearray)
end

#let's also plot this-- we will want to come back here when we do the bilinear-biquadratic model so we know we're getting the same ground state.
plot(2:Ns,params.SvNvec, label = false, ylabel = "von Neumann Entanglement Entropy", xlabel = "Bipartition placement (between site i and i+1)")







###Now let's play around with some correlation functions

id, X, Y, YL, YR, Z, BLL, BLR, BQL, BQR = spinOps()
correlationfuncs = zeros(Float64,Ns,Ns)
k=10
for i in 1:Ns, j in 1:Ns
    cor  = mpoterm(1.,[X,X],[i,j],base)+ mpoterm(1.,[YL,YR],[i,j],base)+ mpoterm(1.,[Z,Z],[i,j],base)
    #one can see that there is a remaining bug in baker's code for self-coupling terms. The correlation function of a site with itself should be two but is zero.
    if i != j
        correlationfuncs[i,j] = expect(psi,cor)
    elseif i == j
        correlationfuncs[i,j] = 2.
    end
end
plot(Vector(1:Ns),correlationfuncs[k,:], legend = false, ylabel = string("Correlation Function with site ",k), xlabel = "Site i")
for i in 1:Ns, j in 1:Ns
    cor  = mpoterm(1.,[X,X],[i,j],base)+ mpoterm(1.,[YL,YR],[i,j],base)+ mpoterm(1.,[Z,Z],[i,j],base)
    #one can see that there is a remaining bug in baker's code for self-coupling terms. The correlation function of a site with itself should be two but is zero.
    if i != j
        correlationfuncs[i,j]  = abs(expect(psi,cor))
    elseif i == j
        correlationfuncs[i,j] = 2.
    end
end
plot(Vector(1:Ns),correlationfuncs[10,:], legend = false, ylabel = string("Magnitude of Correlation Function with site ",k), xlabel = "Site i")
