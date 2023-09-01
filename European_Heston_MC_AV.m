% File: European_Heston_MC_AV.m
%
% Purpose: Antithetic Variate variance reduction Monte Carlo 
%          simulations for pricing European Options under the 
%          Heston model
%
% Algorithm: Nicolas Kuiper and Martin Westberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [type, option_price, std_deviation, elapsed_time] = European_Heston_MC_AV(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,T,R)
% set random number generator seed for reproducibility
rng('default'); 
tic
h = T/Nt;
% generate correlated Brownian motion
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
% pre-allocate memory for price and variance paths
X1 = cell(2,1); 
X1{1} = zeros(Nt,Nsim);
X1{2} = zeros(Nt,Nsim);
% antithetics
X2 = cell(2,1);
X2{1} = zeros(Nt,Nsim);         
X2{2} = zeros(Nt,Nsim);        
% initiate asset price and variance at time 0
X1{1}(1,:) = S0;
X1{2}(1,:) = V0;
% antithetics
X2{1}(1,:) = S0;         
X2{2}(1,:) = V0;         
% simulate asset paths (and antithetic paths)
for i = 1:Nt-1
    % generate asset price paths
    X1{1}(i+1,:) = X1{1}(i,:).*exp((r - 0.5*X1{2}(i,:))*h + ...
        sqrt(X1{2}(i,:)).*dW1{1}(i,:));
    % generate asset volatility paths
    X1{2}(i+1,:) = X1{2}(i,:) + kappa*(theta - X1{2}(i,:))*h + ...
        sigma*sqrt(X1{2}(i,:)).*dW1{2}(i,:);
    % ensure volatility is non-negative
    X1{2}(i+1,:) = max(X1{2}(i+1,:), 0);  
    % generate antithetic price paths
    X2{1}(i+1,:) = X2{1}(i,:).*exp((r - 0.5*X2{2}(i,:))*h + ...
        sqrt(X2{2}(i,:)).*dW2{1}(i,:));
    % generate antithetic volatility paths
    X2{2}(i+1,:) = X2{2}(i,:) + kappa*(theta - X2{2}(i,:))*h + ...
        sigma*sqrt(X2{2}(i,:)).*dW2{2}(i,:);
    % ensure volatility is non-negative
    X2{2}(i+1,:) = max(X2{2}(i+1,:), 0);  
end
% define payoff function for option type
if strcmp(type, 'call')
    payoff = max(X1{1}(end,:) - K, 0);
    antithetic_payoff = max(X2{1}(end,:) - K, 0);
elseif strcmp(type,'put')
    payoff = max(K - X1{1}(end,:), 0);
    antithetic_payoff = max(K - X2{1}(end,:), 0);
end
% calculate option price
payoff_mean = mean(payoff);
price = R*payoff_mean;
payoff_std = std(payoff);
% calculate antithetic option price
antithetic_payoff_mean = mean(antithetic_payoff);
antithetic_price = R*antithetic_payoff_mean;
antithetic_payoff_std = std(antithetic_payoff);
% calculate option prices average and time used
option_price = 0.5*(price + antithetic_price);
std_deviation = 0.5*sqrt(payoff_std^2+antithetic_payoff_std^2);
elapsed_time = toc;
end