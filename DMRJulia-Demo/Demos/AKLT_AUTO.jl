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

id, X, Y, YL, YR, Z, BLL, BLR, BQL, BQR = spinOps()
Ns = 5 #start with 5, 40-50 is good longer chain, should take about 10 seconds to construct hamiltonian
MXM = 10
psi = makePsi0(Ns)

Id = Matrix{Float64}(I, 3, 3)
base = Array{Float64,2}[Id for i = 1:Ns]


function blc(N::Int64, a::Number)
    sitecouple1 = mpoterm(a, [X,X],[1,2],base) + mpoterm(a, [YL,YR],[1,2],base) + mpoterm(a, [Z,Z],[1,2],base)
    for i in 2:(N-1)
        sitecouple1 = sitecouple1 + mpoterm(a, [X,X],[i,i+1],base) + mpoterm(a, [YL,YR],[i,i+1],base) + mpoterm(a, [Z,Z],[i,i+1],base)
    end
    return sitecouple1
end
function bqc(N::Int64,a::Number)
    sitecouple = mpoterm(a,[X^2,X^2],[1,2],base) + mpoterm(a, [X*YL,X*YR],[1,2],base) + mpoterm(a, [X*Z,X*Z],[1,2],base) + mpoterm(a, [YL*X,YR*X],[1,2],base) + mpoterm(a, [YL*YL,YR*YR],[1,2],base) + mpoterm(a, [YL*Z,YR*Z],[1,2],base) + mpoterm(a, [Z*X,Z*X],[1,2],base) + mpoterm(a, [Z*YL,Z*YR],[1,2],base) + mpoterm(a, [Z^2,Z^2],[1,2],base)
    for i in 2:N-1
        sitecouple = sitecouple + mpoterm(a,[X^2,X^2],[i,i+1],base) + mpoterm(a, [X*YL,X*YR],[i,i+1],base) + mpoterm(a, [X*Z,X*Z],[i,i+1],base) + mpoterm(a, [YL*X,YR*X],[i,i+1],base) + mpoterm(a, [YL*YL,YR*YR],[i,i+1],base) + mpoterm(a, [YL*Z,YR*Z],[i,i+1],base) + mpoterm(a, [Z*X,Z*X],[i,i+1],base) + mpoterm(a, [Z*YL,Z*YR],[i,i+1],base) + mpoterm(a, [Z^2,Z^2],[i,i+1],base)
    end
    return sitecouple
end



@time H = blc(Ns,1/2)+bqc(Ns,1/6)

@time params = dmrg(psi,H,maxm=MXM,sweeps=500,cutoff=1E-6)
#spits results in a list we have called "params" here, for example
print("Energy = ",params.energy)
print("Max truncation = ", params.maxtrunc)
#with some modification to Baker's code one can also get entanglement features... More on this shortly.
params.Dvec
