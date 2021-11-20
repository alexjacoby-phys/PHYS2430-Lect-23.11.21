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

"""
    Module: Opt

Operator definitions for spin systems, fermions, and t-J model
"""
module Opt
import LinearAlgebra

  function spinOps()
    id = [1 0 0; 0 1 0; 0 0 1]
    Z = [1 0 0 ; 0 0 0; 0 0 -1]
    X = (1/sqrt(2))*[0 1 0; 1 0 1; 0 1 0 ]
    Y = (1/sqrt(2))*[0 -im 0; im 0 -im ; 0 im 0]
    YL = (1/sqrt(2))*[0 -1 0; 1 0 -1 ; 0 1 0]
    YR = -(1/sqrt(2))*[0 -1 0; 1 0 -1 ; 0 1 0]
    BLL = [X,YL,Z]
    BLR = [X,YR,Z]
    BQL = [X^2,X*YL,X*Z,YL*X,YL^2,YL*Z,Z*X,Z*YL,Z^2]
    BQR = [X^2,X*YR,X*Z,YR*X,YR^2,YR*Z,Z*X,Z*YR,Z^2]
    return id, X, Y, YL, YR, Z, BLL, BLR, BQL, BQR
  end
  export  spinOps
end
