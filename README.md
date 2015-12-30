# Multifractal.jl
This library consists in a collection of different methods to calculate the multifractal variables of time series. The library presents methods for the following approaches:

1 - Chhabra-Jensen method

This approach, presented by Chhabra & Jensen (1989) [1], determines the multifractal spectra directly from the signal without using a Legendre transform.

To run an example of this method, you can do the following:

        using Multifractal

        function main(inputfile::ASCIIString, extensionDq::ASCIIString, extensionFa::ASCIIString, extensionTau::ASCIIString, initialQ::Float64, finalQ::Float64, dq::Float64, Np::Int64, r2dq::Float64, r2fa::Float64, scalesToRemove::Int64)

        #Load the data
            data = readdlm(inputfile,' ');
            x = data[:,1];
            y = data[:,2];

            Multifractal.ChhabraJensen(inputfile, extensionDq, extensionFa, extensionTau, x, y, initialQ, finalQ, dq, Np, r2dq, r2fa, scalesToRemove)

        end

        @time main("series.txt","tdq","tfa","tau",-5.0,5.0,1.0,9,-1.0,-1.0,1);

The first parameter of the function main is the input file with the time series you want to study the Multifractal spectrum. We will detail the other parameters later.

Multifractal.ChhabraJensen function returns as outputs 4 different files:

        series.tdq: 
        series.tfa:
        series.tau:
        summaryDq.dat: 

2 - MFDMA

The MFDMA is an approach based detrended moving average (DMA) for multifractal analyses [2].

3 - MFDFA



[3}


  ```
  julia> using Multifractal

  ```

# References

[1] - Chhabra, A., & Jensen, R. V. (1989). Direct determination of the f (α) singularity spectrum. Physical Review Letters, 62(12), 1327.

[2] - Gu, G. F., & Zhou, W. X. (2010). Detrending moving average algorithm for multifractals. Physical Review E, 82(1), 011136.
