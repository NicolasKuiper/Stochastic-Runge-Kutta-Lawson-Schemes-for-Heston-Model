%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File: Asian_Heston_Midpoint_FSL_CV.m
%
% Purpose: Monte Carlo simulation of the Heston model by a Midpoint Full  
%          Stochastic Lawson scheme with control variance reduction 
%          technique to price Asian Options
%
% Algorithm: Kristian Debrabant, Anne Kværnø, Nicky Gordua Matsson.
%            Runge-Kutta Lawson schemes for stochastic differential 
%            equations. BIT Numerical Matematics 61 (2021), 381-409.
%
% Implementation: Kristian Debrabant, Anne Kværnø, Nicky Gordua Matsson.
%                 Matlab code: Runge-Kutta Lawson schemes for stochastic 
%                 differential equations (2020). 
%                 https://doi.org/10.5281/zenodo.4062482 
%
% Adapted by Nicolas Kuiper and Martin Westberg
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [type, arithmetic_price, geometric_price, arithmetic_std, geometric_std, elapsed_time] = Asian_Heston_Midpoint_FSL_CV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R)
tic
% Prepare input parameters to call Matlab function MidpointFSLVectorized
tspan=[0,T]; 
X0=[S0;V0];
ExpMatrixB=cell(2);
ExpMatrixB{1}=@(p,dW) RotMatExpdW(p,dW);
ExpMatrixB{2}=@(p,dW) RotMatExpdW(p,dW);

% Calculate g1 and g2
g{1}=@(x)[sqrt(x(2,:)).*x(1,:);zeros(1,size(x,2))];
g{2}=@(x)[zeros(1,size(x,2));sigma*sqrt(x(2,:))];

% g{1} = @(x)[sqrt(real(x(2,:))).*x(1,:);sigma*rho*sqrt(real(x(2,:)))];       
% g{2} = @(x)[zeros(1,size(x,2));sigma*sqrt((1-(rho^2))*real(x(2,:)))];

gJac{1}=@(x) getgJac1(x);
gJac{2}=@(x) getgJac2(x);
B=cell(2);
Bexp=cell(2);
B{1}=zeros(2);
B{2}=zeros(2);
Bexp{1} = @(W) ExpMatrixB{1}(0,W);
Bexp{2} = @(W) ExpMatrixB{1}(0,W);

h = T/Nt;
% Generate Brownian Motions
%rng('default'); 
Z1 = randn(Nt,Nsim);
Z2 = randn(Nt,Nsim);            
dW = cell(2,1);
dW{1} = sqrt(h)*Z1;
dW{2} = rho*dW{1} + sqrt(h)*sqrt(1 - rho^2)*Z2;

A=[r,0;0,-kappa];
g0=@(x) getg0(x);
g0Jac=@(x) getgJac0(x);

% Generate asset price at maturity
[~,S,~] = MidpointFSLVectorized(X0,A,g0,B,g,tspan,h,dW,g0Jac,gJac,Bexp);    % returns price at expiration, and price and volatility paths
S = real(S);

% Calculate the price of the Asian option
arithmetic_mean = zeros(1,Nsim);
geometric_mean = zeros(1,Nsim);

for i = 1:Nsim
    arithmetic_mean(i) = mean(S(2:end,i));
    geometric_mean(i) = geomean(S(2:end,i));
end

if strcmp(type,'call')
    arithmetic_payoff = max(arithmetic_mean - K,0);      
    geometric_payoff = max(geometric_mean - K,0);
else
    arithmetic_payoff = max(K - arithmetic_mean,0);
    geometric_payoff = max(K - geometric_mean,0);
end
% set the geometric payoff as control variate
CV = geometric_payoff;
payoff = arithmetic_payoff;
% estimate the control variate coefficient
v = var(payoff);
C = cov(CV, payoff);
b = C(1,2)/v;   
adjusted_payoff = payoff - b*(CV - mean(CV));
arithmetic_price = R*mean(adjusted_payoff);
geometric_price = R*mean(geometric_payoff);

arithmetic_std = std(adjusted_payoff);
geometric_std = std(geometric_payoff);
elapsed_time = toc; 

%% Functions for calling Midpoint FSL
function [erg,inverg]=RotMatExpdW(lambda,dW)
%Calculate Matrix exponentials expm( [0 -lambda;lambda 0]*dW(i)) and their
%inverses and save them in erg(:,:,i) and inverg(:,:,i), respectively.
    erg=zeros(2,2,length(dW));
    inverg=zeros(size(erg));
    temp1=cos(lambda*dW);
    temp2=sin(lambda*dW);
    erg(1,1,:)=temp1;
    erg(2,2,:)=temp1;
    erg(1,2,:)=-temp2;
    erg(2,1,:)=temp2;
    inverg(1,1,:)=temp1;
    inverg(2,2,:)=temp1;
    inverg(1,2,:)=temp2;
    inverg(2,1,:)=-temp2;
end

function result=getg0(x)
    result=ones(size(x));
    result(1,:)=0*result(1,:);
    result(2,:)=result(2,:)*kappa*theta;
end

function result=getgJac0(x)
    nw=size(x,1);
    P=size(x,2);
    result=zeros(nw,nw,P);
end

function result=getgJac1(x)
    nw=size(x,1);
    P=size(x,2);
    result=zeros(nw,nw,P);
    result(1,1,:)=sqrt(x(2,:));
    result(1,2,:)=x(1,:)./(2*sqrt(x(2,:)));
end

function result=getgJac2(x)
    nw=size(x,1);
    P=size(x,2);
    result=zeros(nw,nw,P);
    result(2,2,:)=sigma./(2*sqrt(x(2,:)));
end

end
