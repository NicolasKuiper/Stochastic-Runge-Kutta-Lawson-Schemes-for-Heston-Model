% File: European_Heston_MC_CV.m
%
% Purpose: Control Variate variance reduction Monte Carlo 
%          simulations for pricing European Options under
%          the Heston model
%
% Algorithm: Nicolas Kuiper and Martin Westberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [type, option_price, std_deviation, elapsed_time] = European_Heston_MC_CV(S0,r,V0,K,K_cv,type,kappa,theta,sigma,rho,Nt,Nsim,T,R)
% set random number generator seed for reproducibility
%rng('default');
tic
h = T/Nt;
% generate correlated Brownian motions
% generate two matrices of standard normal numbers
Z1 = randn(Nt,Nsim);           
Z2 = randn(Nt,Nsim);            
dW1 = cell(2,1);
% calculate first Brownian motion matrix
dW1{1} = sqrt(h)*Z1;         
% calculate correlated Brownian motion
dW1{2} = rho*dW1{1} + sqrt(h)*sqrt(1 - rho^2)*Z2; 
% pre-allocate memory for price and variance paths
X1 = cell(2,1); 
X1{1} = zeros(Nt,Nsim);
X1{2} = zeros(Nt,Nsim);
% initiate asset price and variance at time 0
X1{1}(1,:) = S0;
X1{2}(1,:) = V0;
% simulate asset paths
for i = 1:Nt-1
    % generate asset price paths
    X1{1}(i+1,:) = X1{1}(i,:).*exp((r - 0.5*X1{2}(i,:))*h + ...
        sqrt(X1{2}(i,:)).*dW1{1}(i,:));
    % generate asset volatility paths
    X1{2}(i+1,:) = X1{2}(i,:) + kappa*(theta - X1{2}(i,:))*h + ...
        sigma*sqrt(X1{2}(i,:)).*dW1{2}(i,:);
    % ensure volatility is positive
    X1{2}(i+1,:) = max(X1{2}(i+1,:), 0);
end
% define payoff function for option type
if strcmp(type, 'call')
    payoff = max(X1{1}(end,:) - K,0);
    payoff_CV = max(X1{1}(end,:) - K_cv, 0);
elseif strcmp(type,'put')
    payoff = max(K - X1{1}(end,:),0);
    payoff_CV = max(K_cv - X1{1}(end,:), 0);
end
% estimate the control variate coefficient
CV = payoff_CV;
v = var(payoff);
C = cov(payoff, CV);
b = C(1,2)/v;
% calculate option price and time used
adjusted_payoff = payoff - b*(CV - mean(CV));
option_price = R*mean(adjusted_payoff);
std_deviation = std(payoff);
elapsed_time = toc;
end