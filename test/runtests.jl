using Multifractal 

function main(inputfile::ASCIIString, extensionDq::ASCIIString, extensionFa::ASCIIString, initialQ::Float64, finalQ::Float64, dq::Float64, Np::Int32, r2dq::Float64, r2fa::Float64, scalesToRemove::Int32)
    Multifractal.ChhabraJansen(inputfile, extensionDq, extensionFa, initialQ, finalQ, dq, Np, r2dq, r2fa, scalesToRemove)
end

#main(ARGS[1], ARGS[2], ARGS[3], parse(Float64,ARGS[4]), parse(Float64,ARGS[5]), parse(Float64,ARGS[6]), parse(Int,ARGS[7]), parse(Float64,ARGS[8]), parse(Float64,ARGS[9]),parse(Int,ARGS[10]));
main("series.txt","tdq","tfa",-5.0,5.0,1.0,9,-1.0,-1.0,0);
