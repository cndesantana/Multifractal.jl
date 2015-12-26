# version 0.1

module Multifractal 

using Vega

function ChhabraJansen(ASCIIString::inputfile,ASCIIString::extensionDq,ASCIIString::extensionFa,Float64::initialQ,Float64::finalQ,Float64::dq,Int64::Np,Float64::r2dq,Float64::r2fa,Int64::typeofoutput,Int64::scalesToRemove)
    println("Chhabra-Jansen multifractal method!");
end

function fitting(Array{Float64,1}::vx, Array{Float64,1}::vy, Int64::N)
    println("Fitting function");
end

function calcSumM(Array{Float64,1}::x,Array{Float64,1}::y,Float64::Ei,Float64::Ef,Int64::N)
    println("calcSum function");
end

#Write the functions here

end #module
