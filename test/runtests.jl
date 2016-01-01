using Multifractal 
using Benchmarks

function main(inputfile::ASCIIString, extensionDq::ASCIIString, extensionFa::ASCIIString, extensionTau::ASCIIString, initialQ::AbstractFloat, finalQ::AbstractFloat, dq::AbstractFloat, Np::Integer, r2dq::AbstractFloat, r2fa::AbstractFloat, scalesToRemove::Integer)

#;#Load the data
    data = readdlm(inputfile,' ');
    x = data[:,1];
    y = data[:,2];

    Multifractal.ChhabraJensen(inputfile, extensionDq, extensionFa, extensionTau, x, y, initialQ, finalQ, dq, Np, r2dq, r2fa, scalesToRemove)

end

#main(ARGS[1], ARGS[2], ARGS[3], "tau", parse(Float64,ARGS[4]), parse(Float64,ARGS[5]), parse(Float64,ARGS[6]), parse(Int,ARGS[7]), parse(Float64,ARGS[8]), parse(Float64,ARGS[9]),parse(Int,ARGS[10]));

@benchmark main("eeg.txt","tdq","tfa","tau",-5.0,5.0,1.0,9,-1.0,-1.0,1);
