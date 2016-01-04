#Great work, Lucas! :) That is the way to do it!
#I am commenting the function just to test some changes I am doing. 
#

using GLM
using DataFrames
using Vega

data = readdlm("brown.txt",' ');
x = data[:,2];
q_ = collect(-5.0:0.1:5.0);
x,n_min,n_max,N,theta,q = x,10,100,30,0,q_
x,n_min,n_max,N,theta,q = x,10,100,30,0,q_

function float_to_integer(n)
   aux = Array(Integer, 0);
   for(i in n)
      push!(aux,convert(Integer,i));
   end
   return aux; 
end

function MFDMA(x,n_min,n_max,N,theta,q)
    M = length(x);
    MIN = log10(n_min);
    MAX = log10(n_max);

    n = unique(round(logspace(MIN,MAX,N)))
    n = float_to_integer(n); 

    # To build a cumulative sum of the vector y
    y = cumsum(x);
    
    for (i in 1:length(n))
        lgth = n[i];
        # Moving average function 
        y1 = zeros(M-lgth+1);
        for (j in 1:(M-lgth+1))
            y1[j] = mean(y[j:(j+lgth-1)]);
        end #end-forj
#       Determine the residual e
        residuals_=y[max(1,floor(lgth*(1-theta))):(max(1,floor(lgth*(1-theta)))+length(y1)-1)]-y1;
#       Estimate the root-mean-square function F
        F=[];#initialize the variable F
        for (k in 1:convert(Integer,(floor(length(residuals_)/lgth))))
            push!(F,sqrt(mean(residuals_[(k-1)*lgth+1:k*lgth].^2)));#to fulfill a vector dynamically you can use the function push!
        end#end-fork

#Calculate the q-th order overall fluctuation function Fq
    Fq = zeros(length(F),length(q));#initializing the matrix Fq
    for (k in 1:length(q))
        for (j in 1:length(F))
            f=F[j];#is 'f' a scalar value or a vector?
            if q[k] == 0
                Fq[j,k]=exp(0.5*mean(log(f.^2)));#we use .^ when we are working with vectors. But, as far as I understand, f is a value and not a vector. If I am right here, why do we need to calculate the mean of a scalar value?  
            else
                Fq[j,k]=(mean(f.^q[k]))^(1/q[k]);#as far as I understand, both, f and q[i], are scalar values. So I assume we don't need to use .^. Also, I don't get why we should calculate a mean here, if the parameter of the mean function is a scalar value.
            end
        end
    end



#####I AM HERE
#    Calculate the multifractal scaling exponent tau(q)
    for (k in 1:length(q)) #should we use 'length(q)' instead of 'size(Fq,2)'?

# I don't know!

    	fq = Fq[:,k];
    	data = DataFrames.DataFrame(log(fq),log(n));
    	OLS = GLM.glm(Y~X,data,Normal(),IdentityLink());#who are 'Y' and 'X'? Shouldn't be 'y' and 'x' instead? 

#It is a parameter of the function glm. It means that I want a function of the type y = ax + b on my regressive model.
        res = coef(OLS);
      	h[k]=res[2];
    end
    
    tau=h.*q-1;

#    Calculate the singularity strength function alpha(q) and spectrum f(alpha) 
    dx=7;
    dx=fix((dx-1)/2);
    for (i in (dx+1):(length(tau)-dx))
       xx=q[i-dx:i+dx];
       yy=tau[i-dx:i+dx];
    	data = DataFrame(xx,yy);
       OLS = glm(Y~X,data,Normal(),IdentityLink());
       res = coef(OLS);
    	alpha[i,1]=res[2];
    end
    alpha=alpha[dx+1:end];
    f=q[dx+1:end-dx].*alpha - tau[dx+1:end-dx];
end

