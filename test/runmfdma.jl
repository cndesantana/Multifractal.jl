using Multifractal

function main(x,n_min,n_max,N,theta,q)

Multifractal.MFDMA(x,n_min,n_max,N,theta,q)

end

data = readdlm("brown.txt",' ');
x = data[:,2];
q_ = collect(-5.0:0.1:5.0);

main (x,10,100,30,0,q_)
