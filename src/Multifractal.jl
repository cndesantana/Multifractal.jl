# version 0.1

module Multifractal 

using Vega

function ChhabraJansen(inputfile::ASCIIString, extensionDq::ASCIIString, extensionFa::ASCIIString, initialQ::Float64, finalQ::Float64, dq::Float64, Np::Int32, r2dq::Float64, r2fa::Float64, scalesToRemove::Int32)
    println("Chhabra-Jansen multifractal method!");
end

function fitting(vx::Array{Float64,1}, vy::Array{Float64,1}, N::Int32)
    println("Fitting function");
end

function calcSumM(x::Array{Float64,1}, y::Array{Float64,1}, Ei::Float64, Ef::Float64, N::Int32)
    println("CalcSum function");
end

#Write the functions here

end #module
