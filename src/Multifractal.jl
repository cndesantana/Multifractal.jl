# version 0.1

module Multifractal 

using Vega

type Hstc 
    sl::AbstractFloat
    sd::AbstractFloat
    r::AbstractFloat
    in::AbstractFloat
    ea::AbstractFloat
    eb::AbstractFloat
end
Hstc(sl,sd,r,in,ea,eb) = Hstc(sl,sd,r,in,ea,eb);

function printPartitionFunction{T}(FoutTau::IOStream, Qi::T, Qf::T, dq::T, Np::Integer, mye::Vector{T}, Md::Matrix{T})

    line = "Scale";
    @inbounds @simd for(q in Qi:dq:Qf) 
       line = string(line," ",q);
    end
    writedlm(FoutTau,[line],' ');

    @inbounds @simd for(k in 1:Np)
        line = string(mye[k]);
        @inbounds @simd for(q in Qi:dq:Qf) 
            line = string(line," ",Md[round(Int,(q-Qi)/dq)+1,k]);
        end
        writedlm(FoutTau,[line],' ');
    end
end

function Chext(filename,extension)
    return(split(filename,'.')[1]*"."*extension);
end

function MFDFA()
end

#Great work, Lucas! :) That is the way to do it!
#I am commenting the function just to test some changes I am doing. 
#
#function MFDMA(x,n_min,n_max,N,theta,q)
#       M = lenght(x);
#       MIN = log10(n_min);
#       Max = log10(n_max);
## n = (unique(round(logspace(MIN,MAX,N)))' to translate
#
## To build a cumulative sum of the vector y
#
#	y = cumsum(x);
#
#	for (i in 1:length(n))
#
#		lgth = n(i,1);
#
## Moving average function 
#
#		y1 = zeros(1,M-lgth+1);
#		for (j in 1:M-lgth+1)
#			y1 (j) = mean(y(j:j+lgth-1);
#		end
#	end
#
#end

function fitting{T}(vx::Vector{T}, vy::Vector{T}, N::Integer)

    for x = [:s,:sx,:sy,:sx2,:sxy,:sy2,:a,:b,:r,:rx,:ry,:w,:sa,:sb]
        @eval $x = AbstractFloat;
    end
    sx=0.0;
    sy=0.0;
    sx2=0.0;
    sxy=0.0;
    sy2=0.0;

    sx = sum(vx);
    sy = sum(vy);
    sxy = sum(vx.*vy);
    sx2 = sum(vx.*vx);
    sy2 = sum(vy.*vy);

    s = sx2 - sx.*sx/N;
    a = (sxy - sx.*sy/N)/s;
    b = (sy - a*sx)/N;
    w = sy2 + a*a*sx2 + N*b*b;
    w = w - 2.0*a*sxy - 2.0*b*sy + 2.0*a*b*sx;
    if (w < 0.0) 
        w = 0.0;
    else 
        w = sqrt(w/(N-2));
    end
    rx = sx2-sx2/N;
    ry = sy2-sy2/N;

#;#    // Slope error
    sa = (sy2 + N*b*b + a*a*sx2 - 2*(b*sy-a*b*sx+a*sxy))/(N-2);
    sb = sqrt( (sx2*sa)/(N*sx2-sx*sx) );
    sa = sqrt( (N*sa)/(N*sx2-sx*sx) );

    if(abs(ry)<1.0e-10)
        if(abs(a)<1.0e-10) 
            r = 1.0;
        else 
            r = 30000.0;
        end
    else 
        r = a*a*rx/ry;
    end
    return Hstc(a,w,r,b,sa,sb)
end

function calcSumM{T}(x::Vector{T}, y::Vector{T}, Ei::T, Ef::T, N::Integer)
    mysum=0.0::AbstractFloat;
    i=1
    @inbounds while i<=N && x[i]<=Ei i+=1 ; end
    j=i
    @inbounds while j<=N && x[j]<=Ef j+=1 ; end
    @inbounds @simd for k=i:(j-1) mysum += y[k] ; end
    return(mysum);
end 

#function calcSumM2(x::Vector{T}, y::Vector{T}, Ei::T, Ef::T, N::Integer)
#    ret=0.0::AbstractFloat;
#    for(i in 1:N)
#        if( Ei < x[i] <= Ef)
#            ret += y[i];
#        end
#    end
#    return ret;
#end
#   
#function calcSumM3{T}(x::Vector{T}, y::Vector{T}, Ei::T, Ef::T, N::Integer)
#    mysum = zero(T)
#    @inbounds @simd for i in eachindex(x, y)
#         mysum += ifelse(Ei < x[i] <= Ef, y[i], zero(T)) 
#         
#    end
#    return mysum
#end 
 
function getMultifractalCoefficients(FAq::Hstc, FFq::Hstc, FDq::Hstc, q, dq, Dq, RmFa, RmDq, Fout::IOStream, FoutFa::IOStream)

    AlphaMin=QAlphaMin=Fmn=Dqmn=999.9::AbstractFloat;
    AlphaMax=QAlphaMax=Fmx=Dqmx=-999.9;
    qMin=EDqmn=RDqmn=EAlphaMin=RAlphaMin=0.0::AbstractFloat;	#// Alfa minimo, erro e r2
    qMax=EDqmx=RDqmx=EAlphaMax=RAlphaMax=0.0::AbstractFloat;	#// Alfa maximo, erro e r2
    Alpha0=EAlpha0=RAlpha0=D0=RD0=ED0=0.0::AbstractFloat;
    D2=D1=RD1=RD2=ED1=ED2=-1.0;	#// -1 indicates that for the especific q (2 or 1) the R was not calculated

    if((FAq.r >= RmFa) && (FFq.r >= RmFa))
       writedlm(FoutFa,[FAq.sl FAq.sd FAq.r FFq.sl FFq.sd FFq.r],' ');
       if(FAq.sl > AlphaMax) 
           AlphaMax = FAq.sl;
           EAlphaMax = FAq.sd;
           RAlphaMax = FAq.r;
           QAlphaMax = q;
       end 
       if(FAq.sl < AlphaMin) 
           AlphaMin = FAq.sl;
           EAlphaMin = FAq.sd;
           RAlphaMin = FAq.r;
           QAlphaMin = q;
       end 
       if(FFq.sl < Fmn) 
           Fmn = FFq.sl;
       end 
       if(FFq.sl > Fmx) 
           Fmx = FFq.sl;
       end 
       if((0-dq/2) < q <(0+dq/2))
           Alpha0 = FAq.sl;
           EAlpha0 = FAq.sd;
           RAlpha0 = FAq.r;
       end 
    end 
    if(FDq.r >= RmDq)
       writedlm(Fout,[q Dq Dq*(q-1) FDq.sd FDq.r],' ');
       if ((1-dq/2) < q <(1+dq/2))
           EDq = FDq.ea
       else
           EDq = abs(FDq.ea/(q-1));
       end
       if(Dq > Dqmx)                                                                                                
          Dqmx = Dq;
          qMax = q;
          EDqmx = EDq;
          RDqmx = FDq.r;
       end
       if(Dq < Dqmn)           
          Dqmn = Dq;
          qMin = q;
          EDqmn = EDq;
          RDqmn = FDq.r;
       end
       if((0-dq/2) < q < (0+dq/2))
           D0 = Dq;
           RD0 = FDq.r;
           ED0 = EDq;
       end
       if((1-dq/2) < q < (1+dq/2))
           D1 = Dq;
           RD1 = FDq.r;
           ED1 = EDq;
       end
       if((2-dq/2) < q < (2+dq/2))
           D2 = Dq;
           RD2 = FDq.r;
           ED2 = EDq;
       end
    end

    return AlphaMin, AlphaMax, QAlphaMax, QAlphaMin, Fmx, Fmn, Dqmx, Dqmn, qMin, qMax, EDqmx, RDqmx, EDqmn, RDqmn, EAlphaMin, RAlphaMin, EAlphaMax, RAlphaMax, D0, RD0, ED0, D1, RD1, ED1, D2, RD2, ED2, Alpha0, EAlpha0, RAlpha0

end

function ChhabraJensen{T}(inputfile::ASCIIString, extensionDq::ASCIIString, extensionFa::ASCIIString, extensionTau::ASCIIString, x::Vector{T}, y::Vector{T}, Qi::T, Qf::T, dq::T, Np::Integer, RmDq::T, RmFa::T, Io::Integer)
    
    NFout = Chext(inputfile,extensionDq);
    NFoutFA = Chext(inputfile,extensionFa);
    NFoutTau = Chext(inputfile,extensionTau);
    NFoutSumm = "summaryDq.dat";
    Fout = open(NFout,"w+");
    FoutFa = open(NFoutFA,"w+");
    FoutTau = open(NFoutTau,"w+");
    FoutSumm = open(NFoutSumm,"a+");

    AlphaMin=QAlphaMin=Fmn=Dqmn=999.9::AbstractFloat;
    AlphaMax=QAlphaMax=Fmx=Dqmx=-999.9;
    qMin=EDqmn=RDqmn=EAlphaMin=RAlphaMin=0.0::AbstractFloat;	#// Alfa minimo, erro e r2
    qMax=EDqmx=RDqmx=EAlphaMax=RAlphaMax=0.0::AbstractFloat;	#// Alfa maximo, erro e r2
    Alpha0=EAlpha0=RAlpha0=D0=RD0=ED0=0.0::AbstractFloat;
    D2=D1=RD1=RD2=ED1=ED2=-1.0;	#// -1 indicates that for the especific q (2 or 1) the R was not calculated

#;#    /* Fix the size of the file, the maximum and minimum */
    Md = zeros(round(Int,((Qf-Qi)/dq)+1),Np+1);
    Ma = zeros(Np+1);
    Mf = zeros(Np+1);
    mye = zeros(Np+1);

    N = length(y);
    MaxY = maximum(y);
    MaxX = maximum(x);
    MinY = minimum(y);
    MinX = minimum(x);
    SomaY = sum(y);
    x = (x-MinX)/(MaxX-MinX);
#;############################# To change from here on
#;    // Begins the "thing"
    I=Io;			#// Initial partition, for I=1 the mi(Epson) finalize with Epson=1/2

        @inbounds @simd for(q in Qi:dq:Qf)
        Md = zeros(round(Int,((Qf-Qi)/dq)+1),Np+1);
        Ma = zeros(Np+1);
        Mf = zeros(Np+1);
        @inbounds @simd for(k in I:Np)						#// Loop for partition numbers
            Nor=0.0::AbstractFloat;
            m=0.0::AbstractFloat;
            Pr=0::Integer;
            Pr = 2^(k-1);
            E = 1.0/Pr;						#// Size of each partition
            mye[k-I+1] = log10(E);
            pos = k-I+1;
            val = mye[pos];

            @inbounds @simd for(i in 1:Pr)						#// To estimate f(alfa)
                m = calcSumM(x,y,(i-1)*E,i*E,N)/SomaY;
                if(m!=0)
                    Nor += m^q;
                end
            end
            
            @inbounds @simd for(i in 1:Pr) #// loop for scan over the partition
                m = calcSumM(x,y,(i-1)*E,i*E,N)/SomaY;
                if(m!=0)		        #// Evita divergencias de medidas nulas
                    currentval = Md[round(Int,(q-Qi)/dq)+1,k-I+1]::AbstractFloat;
                    if( (1-dq/2) < q < (1+dq/2) )
                        setindex!(Md,currentval + m*log10(m)/Nor,round(Int,(q-Qi)/dq)+1,k-I+1)
                    else    
                        setindex!(Md,currentval + m^q,round(Int,(q-Qi)/dq)+1,k-I+1)
                    end
                    mq = (m^q)/Nor;					#// To estimate f(alfa)
                    currentval = Ma[k-I+1]::AbstractFloat;
                    setindex!(Ma,currentval + mq*log10(m),k-I+1);
                    currentval = Mf[k-I+1]::AbstractFloat;
                    setindex!(Mf,currentval + mq*log10(mq),k-I+1);
                end #end-if
            end #end-for

            if(! ((1-dq/2) < q < (1+dq/2)) )
                setindex!(Md,log10(Md[round(Int,(q-Qi)/dq)+1,k-I+1]),round(Int,(q-Qi)/dq)+1,k-I+1); #// if q!=1
            end        
        end        

        FAq = fitting(mye,Ma,Np);
        FFq = fitting(mye,Mf,Np);
        FDq = fitting(mye,Md'[:,round(Int,(q-Qi)/dq)+1],Np);
        if( (1-dq/2) < q < (1+dq/2) )
            Dq = FDq.sl::AbstractFloat;
        else 
            Dq = FDq.sl/(q-1)::AbstractFloat;
            FDq.sd /= abs(q-1);
        end

#        AlphaMin, AlphaMax, QAlphaMax, QAlphaMin, Fmx, Fmn, Dqmx, Dqmn, qMin, qMax, EDqmx, RDqmx, EDqmn, RDqmn, EAlphaMin, RAlphaMin, EAlphaMax, RAlphaMax, D0, RD0, ED0, D1, RD1, ED1, D2, RD2, ED2, Alpha0, EAlpha0, RAlpha0 = getMultifractalCoefficients(FAq, FFq, FDq, q, dq, Dq, RmFa, RmDq, Fout, FoutFa);

        if((FAq.r >= RmFa) && (FFq.r >= RmFa))
           writedlm(FoutFa,[FAq.sl FAq.sd FAq.r FFq.sl FFq.sd FFq.r],' ');
           if(FAq.sl > AlphaMax) 
               AlphaMax = FAq.sl;
               EAlphaMax = FAq.sd;
               RAlphaMax = FAq.r;
               QAlphaMax = q;
           end 
           if(FAq.sl < AlphaMin) 
               AlphaMin = FAq.sl;
               EAlphaMin = FAq.sd;
               RAlphaMin = FAq.r;
               QAlphaMin = q;
           end 
           if(FFq.sl < Fmn) 
               Fmn = FFq.sl;
           end 
           if(FFq.sl > Fmx) 
               Fmx = FFq.sl;
           end 
           if((0-dq/2) < q <(0+dq/2))
               Alpha0 = FAq.sl;
               EAlpha0 = FAq.sd;
               RAlpha0 = FAq.r;
           end 
        end 
        if(FDq.r >= RmDq)
           writedlm(Fout,[q Dq Dq*(q-1) FDq.sd FDq.r],' ');
           if ((1-dq/2) < q <(1+dq/2))
               EDq = FDq.ea
           else
               EDq = abs(FDq.ea/(q-1));
           end
           if(Dq > Dqmx)                                                                                                
              Dqmx = Dq;
              qMax = q;
              EDqmx = EDq;
              RDqmx = FDq.r;
           end
           if(Dq < Dqmn)           
              Dqmn = Dq;
              qMin = q;
              EDqmn = EDq;
              RDqmn = FDq.r;
           end
           if((0-dq/2) < q < (0+dq/2))
               D0 = Dq;
               RD0 = FDq.r;
               ED0 = EDq;
           end
           if((1-dq/2) < q < (1+dq/2))
               D1 = Dq;
               RD1 = FDq.r;
               ED1 = EDq;
           end
           if((2-dq/2) < q < (2+dq/2))
               D2 = Dq;
               RD2 = FDq.r;
               ED2 = EDq;
           end
        end
    end

    writedlm(FoutSumm,[inputfile qMin qMax Dqmn EDqmn RDqmn Dqmx EDqmx RDqmx D0 ED0 RD0 D1 ED1 RD1 D2 ED2 RD2 QAlphaMin QAlphaMax Alpha0 EAlpha0 RAlpha0 AlphaMax EAlphaMax RAlphaMax AlphaMin EAlphaMin RAlphaMin Fmn Fmx],'\t');

    printPartitionFunction(FoutTau, Qi, Qf, dq, Np, mye, Md);
    close(Fout);
    close(FoutFa);
    close(FoutTau);
    close(FoutSumm);
end

#Write the functions here

end #module
