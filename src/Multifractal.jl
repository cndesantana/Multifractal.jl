# version 0.1

module Multifractal 

using Vega

function start()
    println("Hello multifractal world!");
end

function fitting(Array{Float64,1}::vx, Array{Float64,1}::vy, Int64::N)
    println("Fitting function");
end

function calc_SumM(Array{Float64,1}::x,Array{Float64,1}::y,Float64::Ei,Float64::Ef,Int64::N)
    println("calcSum function");
end

#Write the functions here

end #module
