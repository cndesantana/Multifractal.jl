% Algorithm developed by Gu, G. F., & Zhou, W. X.

function [n,Fq,tau,alpha,f]=OneDimensional_MFDMA(x,n_min,n_max,N,theta,q)
% The procedure [n,Fq,tau,alpha,f]=OneDimensional_MFDMA(x,n_min,n_max,N,,theta,q) 
% is used to calculate the multifractal properties of one-dimensional time
% series.

% Input:
%   x: the time series we considered
%   n_min: the lower bound of the segment size n
%   n_max: the upper bound of the segment size n
%   N: the length of n, that is, the data points in the plot of Fq VS n
%   theta: the position parameter of the sliding window
%   q: multifractal order 

% Output:
%   n: segment size series
%   Fq: q-th order fluctuation function
%   tau: multifractal scaling exponent
%   alpha: multifractal singularity strength function
%   f: multifractal spectrum

% The procedure works as follows:
%   1) Construct the cumulative sum y.
% 
%   2) For each n, calculate the moving average function \widetilde{y}. 
% 
%   3) Determine the residual e by detrend \widetilde{y} from y.
% 
%   4) Estimate the root-mean-square function F.
% 
%   5) Calculate the q-th order overall fluctuation function Fq.
%
%   6) Calculate the multifractal scaling exponent tau(q).
%
%   7) Calculate the singularity strength function alpha(q) and spectrum
%   f(alpha).

% Note:
%   1) The window size and the segment size must be identical. 
% 
%   2) The lower bound n_min would better be selected around 10.
% 
%   3) The upper bound n_max would better be selected around 
%      ten percent of the length of time series.
% 
%   4) N would better be seleceted in the range [20,40].
%
%   5) The parameter theta varies in the range [0,1]. We recommend theta=0.

% Example: 
%   [n,Fq,tau,alpha,f]=OneDimensional_MFDMA(x,10,round(length(x)/20),30,0,-5:0.1:5);
% Input:
%   x: the time series we considered
%   n_min: the lower bound of the segment size n
%   n_max: the upper bound of the segment size n
%   N: the length of n, that is, the data points in the plot of Fq VS n
%   theta: the position parameter of the sliding window
%   q: multifractal order 


if size(x,2) == 1
    x=x';
end

M=length(x);
MIN=log10(n_min);
MAX=log10(n_max);
n=(unique(round(logspace(MIN,MAX,N))))';


% Construct the cumulative sum y
y=cumsum(x);

for i=1:length(n)
    lgth=n(i,1);
    
    % Calculate the moving average function \widetilde{y}
    y1=zeros(1,M-lgth+1);
    for j=1:M-lgth+1
        y1(j)=mean(y(j:j+lgth-1));
    end
    
    % Determine the residual e
    e=y(max(1,floor(lgth*(1-theta))):max(1,floor(lgth*(1-theta)))+length(y1)-1)-y1;
    
    % Estimate the root-mean-square function F
    for k=1:floor(length(e)/lgth)
        F{i}(k)=sqrt(mean(e((k-1)*lgth+1:k*lgth).^2));
    end
end


% Calculate the q-th order overall fluctuation function Fq
for i=1:length(q)
    for j=1:length(F)
        f=F{j};
        if q(i) == 0
            Fq(j,i)=exp(0.5*mean(log(f.^2)));
        else
            Fq(j,i)=(mean(f.^q(i)))^(1/q(i));
        end
    end
end


% Calculate the multifractal scaling exponent tau(q)
for i=1:size(Fq,2)
    fq=Fq(:,i);
    r=regstats(log(fq),log(n),'linear',{'tstat'});
    k=r.tstat.beta(2);
    h(i,1)=k;
end
tau=h.*q'-1;


% Calculate the singularity strength function alpha(q) and spectrum f(alpha) 
dx=7;
dx=fix((dx-1)/2);
for i=dx+1:length(tau)-dx
    xx=q(i-dx:i+dx);
    yy=tau(i-dx:i+dx);
    r=regstats(yy,xx,'linear',{'tstat'});
    alpha(i,1)=r.tstat.beta(2);
end
alpha=alpha(dx+1:end);
f=q(dx+1:end-dx)'.*alpha-tau(dx+1:end-dx);

%save ('res.dat','n','Fq','tau','alpha','f','F(q)', 'f(alfa)', '-ascii')
save ('res_n.dat','n', '-ascii')
save ('res_Fq.dat','Fq', '-ascii')
save ('res_tau.dat','tau', '-ascii')
save ('res_alpha.dat','alpha', '-ascii')
save ('res_f.dat','f', '-ascii')

figure
plot(q,Fq)
xlabel('q')
ylabel('Fq')

figure
plot(alpha,f)
xlabel('alpha')
ylabel('f')
