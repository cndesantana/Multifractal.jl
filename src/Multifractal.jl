# version 0.1

module Multifractal 

using Vega

type Hstc 
    sl::Float64
    sd::Float64
    r::Float64
    in::Float64
    ea::Float64
    eb::Float64
end
Hstc(sl,sd,r,in,ea,eb) = Hstc(sl,sd,r,in,ea,eb);

function Chext(filename,extension)
    return(split(filename,'.')[1]*"."*extension);
end

function MFDFA()
end

function MFDMA()
end

function fitting(vx::Array{Float64,1}, vy::Array{Float64,1}, N::Int64)
    println("Fitting function");

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

function calcSumM(x::Array{Float64,1}, y::Array{Float64,1}, Ei::Float64, Ef::Float64, N::Int64)
    ret::Float64 = 0;
    for(i in 1:N)
        if( Ei < x[i] <= Ef)
            ret += y[i];
        end
    end
    return ret;
end
    
function getMaxMinFitting(FAq::Hstc, FFq::Hstc, FDq::Hstc, dq::Float64, q::Float64, Dq::Float64, RmFa::Float64, RmDq::Float64, Fout::IOStream, FoutFa::IOStream)
    Amn=999;  
    Amx=-999; 
    Aqmx=-999;
    Aqmn=999;
    Fmx=-999; 
    Fmn=999;  
    Dqmx=-999.9;
    Dqmn= 999.9;

    qMin::Float64;
    qMax::Float64;
    EDqmn::Float64;
    RDqmn::Float64;
    EDqmx::Float64;
    RDqmx::Float64;
    EAmn::Float64;
    RAmn::Float64;	#// Alfa minimo, erro e r2
    EAmx::Float64;
    RAmx::Float64;	#// Alfa maximo, erro e r2
    myDo::Float64;
    RDo::Float64;
    EDo::Float64;

    if((FAq.r >= RmFa) && (FFq.r >= RmFa))
       writedlm(FoutFa,[FAq.sl FAq.sd FAq.r FFq.sl FFq.sd FFq.r],' ');
       if(FAq.sl > Amx) 
           Amx = FAq.sl;
           EAmx = FAq.sd;
           RAmx = FAq.r;
           Aqmx = q;
       end 
       if(FAq.sl < Amn) 
           Amn = FAq.sl;
           EAmn = FAq.sd;
           RAmn = FAq.r;
           Aqmn = q;
       end 
       if(FFq.sl < Fmn) 
           Fmn = FFq.sl;
       end 
       if(FFq.sl > Fmx) 
           Fmx = FFq.sl;
       end 
       if((0-dq/2) < q <(0+dq/2))
           ao = FAq.sl;
           EAo = FAq.sd;
           RAo = FAq.r;
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
           myDo = Dq;
           RDo = FDq.r;
           EDo = EDq;
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
    return Amn, Amx, Aqmx, Aqmn, Fmx, Fmn, Dqmx, Dqmn, qMin, qMax, EDqmx, RDqmx, EDqmn, RDqmn, EAmn, RAmn, EAmx, RAmx, myDo, RDo, EDo

end

function ChhabraJensen(inputfile::ASCIIString, extensionDq::ASCIIString, extensionFa::ASCIIString, Qi::Float64, Qf::Float64, dq::Float64, Np::Int64, RmDq::Float64, RmFa::Float64, Io::Int64)
    
    D2=D1=RD1=RD2=ED1=ED2=-1;	#// -1 indicates that for the especific q (2 or 1) the R was not calculated

    NFout = Chext(inputfile,extensionDq);
    NFoutFA = Chext(inputfile,extensionFa);
    Fout = open(NFout,"w+");
    FoutFa = open(NFoutFA,"w+");

#;#    /* Fix the size of the file, the maximum and minimum */
    Md = zeros(round(Int,((Qf-Qi)/dq)+1),Np+1);
    Ma = zeros(Np+1);
    Mf = zeros(Np+1);
    e = zeros(Np+1);

#;#    /* Load data  */
    data = readdlm(inputfile,' ');
    N = length(data[:,1]);
    x = data[:,1];
    y = data[:,2];
    MaxY = maximum(y);
    MaxX = maximum(x);
    MinY = minimum(y);
    MinX = minimum(x);
    SomaY = sum(y);
    x = (x-MinX)/(MaxX-MinX);
#;############################# To change from here on
#;    // Begins the "thing"
    I=Io;			#// Initial partition, for I=1 the mi(Epson) finalize with Epson=1/2
    for(q in Qi:dq:Qf)
        for(k in I:Np)						#// Loop for partition numbers
            println("k = $k");
            Nor=0.0::Float64;
            m=0.0::Float64;
            Pr=0::Int64;
            Pr = 2^k;
            E = 1.0/Pr;						#// Size of each partition
            e[k-I+1] = log10(E);

            for(i in 1:Pr)						#// To estimate f(alfa)
                println("i1 = $i");
                m = calcSumM(x,y,(i-1)*E,i*E,N)/SomaY;
                if(m!=0)
                    Nor += m^q;
                end
            end
            
            for(i in 1:Pr) #// loop for scan over the partition
                println("i2 = $i");
                m = calcSumM(x,y,(i-1)*E,i*E,N)/SomaY;
                if(m!=0)		        #// Evita divergencias de medidas nulas
                    currentval = Md[round(Int,(q-Qi)/dq)+1,k-I+1]::Float64;
                    if( (1-dq/2) < q < (1+dq/2) )
                        setindex!(Md,currentval + m*log10(m)/Nor,round(Int,(q-Qi)/dq)+1,k-I+1)
                    else    
                        setindex!(Md,currentval + m^q,round(Int,(q-Qi)/dq)+1,k-I+1)
                    end
                    mq = (m^q)/Nor;					#// To estimate f(alfa)
                    currentval = Ma[k-I+1]::Float64;
                    setindex!(Ma,currentval + mq*log10(m),k-I+1);
                    currentval = Mf[k-I+1]::Float64;
                    setindex!(Mf,currentval + mq*log10(mq),k-I+1);
                end
            end

            if(! ((1-dq/2) < q < (1+dq/2)) )
                setindex!(Md,log10(Md[round(Int,(q-Qi)/dq)+1,k-I+1]),round(Int,(q-Qi)/dq)+1,k-I+1); #// if q!=1
            end        
        end        

        FAq = fitting(e,Ma,Np);
        FFq = fitting(e,Mf,Np);
        FDq = fitting(e,Md'[:,round(Int,(q-Qi)/dq)+1],Np);
        if( (1-dq/2) < q < (1+dq/2) )
            Dq = FDq.sl::Float64;
        else 
            Dq = FDq.sl/(q-1)::Float64;
            FDq.sd /= abs(q-1);
        end

        Amn, Amx, Aqmx, Aqmn, Fmx, Fmn, Dqmx, Dqmn, qMin, qMax, EDqmx, RDqmx, EDqmn, RDqmn, EAmn, RAmn, EAmx, RAmx, myDo, RDo, EDo = getMaxMinFitting(FAq, FFq, FDq, dq, q, Dq, RmFa, RmDq, Fout, FoutFa);
        writedlm(STDOUT,[inputfile qMin qMax Dqmn EDqmn RDqmn Dqmx EDqmx RDqmx myDo EDo RDo D1 ED1 RD1 D2 ED2 RD2 Aqmn Aqmx ao EAo RAo Amx EAmx RAmx Amn EAmn RAmn Fmn Fmx],'\t');
    end

    println("Chhabra-Jansen multifractal method!");
end

#Write the functions here

end #module
