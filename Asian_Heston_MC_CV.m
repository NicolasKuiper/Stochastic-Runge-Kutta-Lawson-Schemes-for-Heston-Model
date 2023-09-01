% File: Asian_Heston_MC_CV.m
%
% Purpose: Monte Carlo simulations with control variate variance 
%          reduction technique for pricing Asian Options under the 
%          Heston model
%
% Algorithm: Nicolas Kuiper and Martin Westberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [type, arithmetic_price,geometric_price, arithmetic_std, geometric_std, elapsed_time] = Asian_Heston_MC_CV(S0, r, V0, K, type, kappa, theta, sigma, rho, Nt, Nsim, T, R)
% set random number generator seed for reproducibility
%rng('default');
tic
h = T/Nt;
% generate correlated Brownian Motion
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
% Initiate asset price and variance at time 0
X1{1}(1,:) = S0;
X1{2}(1,:) = V0;
% simulate asset paths
for i = 1:Nt-1
    % generate asset price paths
    X1{1}(i+1,:) = X1{1}(i,:).*exp((r - 0.5*X1{2}(i,:))*h + sqrt(X1{2}(i,:)).*dW1{1}(i,:));
    % generate asset volatility paths
    X1{2}(i+1,:) = X1{2}(i,:) + kappa*(theta - X1{2}(i,:))*h + sigma*sqrt(X1{2}(i,:)).*dW1{2}(i,:);
    % ensure volatility is non-negative
    X1{2}(i+1,:) = max(X1{2}(i+1,:), 0);
end
% allocate memory for averages
arithmetic_mean = zeros(1,Nsim);
geometric_mean = zeros(1,Nsim);
% calculate the averages for each path
for i = 1:Nsim
    arithmetic_mean(i) = mean(X1{1}(2:end,i));
    geometric_mean(i) = geomean(X1{1}(2:end,i));
end
% define payoff functions and
if strcmp(type,'call')
    arithmetic_payoff = max(arithmetic_mean - K ,0);      
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
% calculate option price
adjusted_payoff = payoff - b*(CV - mean(CV));
arithmetic_price = R*mean(adjusted_payoff);
geometric_price = R*mean(geometric_payoff);
arithmetic_std = std(arithmetic_payoff);
geometric_std = std(geometric_payoff);
elapsed_time = toc;
end