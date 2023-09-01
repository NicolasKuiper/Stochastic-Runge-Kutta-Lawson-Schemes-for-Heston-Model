% File: European_Heston_Midpoint_FSL_CV.m
%
% Purpose: Monte Carlo simulation of the Heston model by a 
%          Midpoint Full Stochastic Lawson scheme with con-
%          trol variate variance reduction technique to price 
%          European Options
%
% Algorithm: Kristian Debrabant, Anne Kv{\ae}rn{\o}, Nicky 
%            Gordua Matsson. Runge-Kutta Lawson schemes for 
%            stochastic differential equations. BIT Numerical 
%            Matematics 61 (2021), 381-409.
%
% Implementation: Kristian Debrabant, Anne Kv{\ae}rn{\o}{\o}, 
%                 Nicky Gordua Matsson.
% Matlab code: Runge-Kutta Lawson schemes for stochastic 
%                 differential equations (2020). 
%                 https://doi.org/10.5281/zenodo.4062482 
%
% Adapted by Nicolas Kuiper and Martin Westberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [type, option_price, std_deviation, elapsed_time] = European_Heston_Midpoint_FSL_CV(S0,r,V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,R)
% set random number generator seed for reproducibility
rng('default');
tic 
h = T/Nt;
% prepare input parameters to call Matlab funxction MidpointFSLVectorized
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
% g1, g2 Jacobians
gJac{1}=@(x) getgJac1(x);
gJac{2}=@(x) getgJac2(x);
% matrices A1, A2
B=cell(2);
Bexp=cell(2);
B{1}=zeros(2);
B{2}=zeros(2);
Bexp{1} = @(W) ExpMatrixB{1}(0,W);
Bexp{2} = @(W) ExpMatrixB{1}(0,W);
% generate Brownian Motions
rng('default'); 
Z1 = randn(Nt,Nsim);
Z2 = randn(Nt,Nsim);            
dW = cell(2,1);
dW{1} = sqrt(h)*Z1;
dW{2} = rho*dW{1} + sqrt(h)*sqrt(1 - rho^2)*Z2;

A=[r,0;0,-kappa];
g0=@(x) getg0(x);
g0Jac=@(x) getgJac0(x);

% Generate asset price at maturity
[~,~,X] = MidpointFSLVectorized(X0,A,g0,B,g,tspan,h,dW,g0Jac,gJac,Bexp);    % returns price at expiration, and price and volatility paths
X = real(X);

% Calculate the price of the European option
if strcmp(type,'call')
    payoff = max(X(1,:) - K,0);
    payoff_CV = max(X(1,:) - K_cv,0);
elseif strcmp(type,'put')
    payoff = max(K - X(1,:),0);
    payoff_CV = max(K_cv - X(1,:),0);
end
CV = payoff_CV;
% Estimate the control variate coefficient
v = var(payoff);
C = cov(payoff, CV);
b = C(1,2)/v;   
adjusted_payoff = payoff - b*(CV - mean(CV));
option_price = R*mean(adjusted_payoff);
std_deviation = std(payoff);
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
