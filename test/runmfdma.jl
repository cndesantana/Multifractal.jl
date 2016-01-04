using Multifractal

function main(x::Array{Float64,1},n_min::Int32,n_max::Int32,
N::Int32,theta::FloatRange{Float64}
,q::Int32)

Multifractal.MFDMA(x,n_min,n_max,N,theta,q)

end

data = readdlm("brown.txt",' ');
x = data[:,2];
q_ = collect(-5.0:0.1:5.0);

main (x,10,100,30,0,q_)
