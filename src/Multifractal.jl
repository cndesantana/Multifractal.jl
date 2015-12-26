# version 0.1

module Multifractal 

using Vega

function start()
    println("Hello multifractal world!");
end

function fitting(vx::Array{Float64,1}, vy::Array{Float64,1}, N::Int64)
    println("Fitting function");
end

function calc_SumM(x::Array{Float64,1},y::Array{Float64,1},Ei::Float64,Ef::Float64,N::Int64)
    println("calcSum function");
end

#Write the functions here

end #module
