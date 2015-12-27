# version 0.1

module Multifractal 

using Vega

function ChhabraJansen(inputfile::ASCIIString, extensionDq::ASCIIString, extensionFa::ASCIIString, initialQ::Float32, finalQ::Float32, dq::Float32, Np::Int32, r2dq::Float32, r2fa::Float32, scalesToRemove::Int32)
    println("Chhabra-Jansen multifractal method!");
end

function fitting(vx::Array{Float32,1}, vy::Array{Float32,1}, N::Int32)
    println("Fitting function");
end

function calcSumM(x::Array{Float32,1}, y::Array{Float32,1}, Ei::Float32, Ef::Float32, N::Int32)
    println("CalcSum function");
end

#Write the functions here

end #module
