% File: Asian_Heston_MC_AV.m
%
% Purpose: Monte Carlo simulations with antithetic variate variance 
%          reduction technique for pricing Asian Options under the Heston 
%          model
%
% Algorithm: Nicolas Kuiper and Martin Westberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [type, arithmetic_price,geometric_price, arithmetic_std, geometric_std, elapsed_time] = Asian_Heston_MC_AV(S0, r, V0, K, type, kappa, theta, sigma, rho, Nt, Nsim, T, R)
% set random number generator seed for reproducibility
rng('default'); 
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
% generate Brownian motions for antithetics path
dW2 = cell(2,1);
dW2{1} = -dW1{1};
dW2{2} = rho*dW2{1} + sqrt(h)*sqrt(1 - rho^2)*(-Z2);
% pre-allocate memory for paths
X1 = cell(2,1); 
X1{1} = zeros(Nt,Nsim);
X1{2} = zeros(Nt,Nsim);
% pre-allocate memory for price and variance paths
X2 = cell(2,1);
X2{1} = zeros(Nt,Nsim);         % antithetics price
X2{2} = zeros(Nt,Nsim);         % antithetics variance
% Initiate asset price and variance at time
X1{1}(1,:) = S0;
X1{2}(1,:) = V0;
X2{1}(1,:) = S0;         % antithetics initial price
X2{2}(1,:) = V0;         % antithetics initial variance
% simulate asset paths under geometric Brownian Motion
for i = 1:Nt-1
    % generate asset price paths
    X1{1}(i+1,:) = X1{1}(i,:).*exp((r - 0.5*X1{2}(i,:))*h + ...
        sqrt(X1{2}(i,:)).*dW1{1}(i,:));
    % generate asset volatility paths 
    X1{2}(i+1,:) = X1{2}(i,:) + kappa*(theta - X1{2}(i,:))*h + ...
        sigma*sqrt(X1{2}(i,:)).*dW1{2}(i,:);
    % ensure volatility is non-negative    
    X1{2}(i+1,:) = max(X1{2}(i+1,:), 0);
    % generate asset price paths for antithetics variate
    X2{1}(i+1,:) = X2{1}(i,:).*exp((r - 0.5*X2{2}(i,:))*h + ...
        sqrt(X2{2}(i,:)).*dW2{1}(i,:));
    % generate asset volatility paths for antithetics variate
    X2{2}(i+1,:) = X2{2}(i,:) + kappa*(theta - X2{2}(i,:))*h + ...
        sigma*sqrt(X2{2}(i,:)).*dW2{2}(i,:);
    % ensure volatility is non-negative
    X2{2}(i+1,:) = max(X2{2}(i+1,:), 0);  
end
% pre-allocate memory for averages
arithmetic_mean = zeros(1,Nsim);
geometric_mean = zeros(1,Nsim);
antithetic_arithmetic_mean = zeros(1,Nsim);
antithetic_geometric_mean = zeros(1,Nsim);
% calculate the averages for each path
for i = 1:Nsim
    arithmetic_mean(i) = mean(X1{1}(2:end,i));
    geometric_mean(i) = geomean(X1{1}(2:end,i));
    antithetic_arithmetic_mean(i) = mean(X2{1}(2:end,i));
    antithetic_geometric_mean(i) = geomean(X2{1}(2:end,i));
end
% calculate payoffs for each path
if strcmp(type,'call')
    arithmetic_payoff = max(arithmetic_mean - K,0);
    geometric_payoff = max(geometric_mean - K,0);
    antithetic_arithmetic_payoff = max(antithetic_arithmetic_mean - K,0);
    antithetic_geometric_payoff = max(antithetic_geometric_mean - K,0);
else
    arithmetic_payoff = max(K - arithmetic_mean,0);
    geometric_payoff = max(K - geometric_mean,0);
    antithetic_arithmetic_payoff = max(K - antithetic_arithmetic_mean,0);
    antithetic_geometric_payoff = max(K - antithetic_geometric_mean,0);
end
% calculate option prices
arithmetic_price = R*mean(arithmetic_payoff);
geometric_price = R*mean(geometric_payoff);
antithetic_arithmetic_price = R*mean(antithetic_arithmetic_payoff);
antithetic_geometric_price = R*mean(antithetic_geometric_payoff);
% average the option prices
arithmetic_price = 0.5*(arithmetic_price + antithetic_arithmetic_price);
geometric_price = 0.5*(geometric_price + antithetic_geometric_price);
arithmetic_std = 0.5*sqrt(std(arithmetic_payoff)^2 + std(antithetic_arithmetic_payoff)^2);
geometric_std = 0.5*sqrt(std(geometric_payoff)^2 + std(antithetic_geometric_payoff)^2); 
elapsed_time = toc;
end