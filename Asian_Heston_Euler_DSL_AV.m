%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File: European_Heston_Euler_DSL_AV.m
%
% Purpose: Monte Carlo simulation of the Heston model by a Euler-
%          Maruyama Drift Stochastic Lawson scheme with antithetic 
%          variate variance reduction technique to price Asian 
%          Options
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
function [type, arithmetic_price, geometric_price, arithmetic_std, geometric_std, elapsed_time] = Asian_Heston_Euler_DSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R)
tic 
X0 = [S0; V0];
A = [r, 0; 0, -kappa];
g0 = @(x) getg0(x,kappa,theta);

g = cell(2, 1);
g{1} = @(x) [sqrt(x(2, :)).*x(1, :); zeros(1, size(x, 2))];
g{2} = @(x) [zeros(1, size(x, 2)); sigma * sqrt(x(2, :))];

tspan = [0, T];
h = T / Nt;
rng('default'); 
Z1 = randn(Nt, Nsim);
Z2 = randn(Nt, Nsim);
dW1 = cell(2, 1);
dW1{1} = sqrt(h) * Z1;
dW1{2} = rho * dW1{1} + sqrt(h) * sqrt(1 - rho^2) * Z2;

% Generate Brownian Motions for antithetic paths
dW2 = cell(2,1);
dW2{1} = -dW1{1};
dW2{2} = rho*dW2{1} + sqrt(h)*sqrt(1 - rho^2)*-Z2;

[~,S1, ~] = EulerDSLVectorized(X0, A, g0, g, tspan, h, dW1);
[~,S2, ~] = EulerDSLVectorized(X0, A, g0, g, tspan, h, dW2);

% Calculate average price throughout option life
arithmetic_mean = zeros(1,Nsim);
antithetic_arithmetic_mean = zeros(1,Nsim);
geometric_mean = zeros(1,Nsim);
antithetic_geometric_mean = zeros(1,Nsim);
for i = 1:Nsim
    arithmetic_mean(i) = mean(S1(2:end,i));
    antithetic_arithmetic_mean(i) = mean(S2(2:end,i));
    
    geometric_mean(i) = geomean(S1(2:end,i));
    antithetic_geometric_mean(i) = geomean(S2(2:end,i));
end

% calculate payoffs for each path
if strcmp(type,'call')
    arithmetic_payoff = max(arithmetic_mean - K,0);
    antithetic_arithmetic_payoff = max(antithetic_arithmetic_mean - K,0);
    
    geometric_payoff = max(geometric_mean - K,0);
    antithetic_geometric_payoff = max(antithetic_geometric_mean - K,0);
else
    arithmetic_payoff = max(k - arithmetic_mean,0);
    antithetic_arithmetic_payoff = max(K - antithetic_arithmetic_mean,0);
    
    geometric_payoff = max(K - geometric_mean,0);
    antithetic_geometric_payoff = max(k - antithetic_geometric_mean,0);
end
% calculate payoffs mean and discount
arithmetic_price = R*mean(arithmetic_payoff);
geometric_price = R*mean(geometric_payoff);

anithetic_arithmetic_price = R*mean(antithetic_arithmetic_payoff);
anithetic_geometric_price = R*mean(antithetic_geometric_payoff);

arithmetic_price = 0.5*(arithmetic_price + anithetic_arithmetic_price);
geometric_price = 0.5*(geometric_price + anithetic_geometric_price);

arithmetic_std = 0.5*sqrt(std(arithmetic_payoff)^2+std(antithetic_arithmetic_payoff)^2);
geometric_std = 0.5*sqrt(std(geometric_payoff)^2+std(antithetic_geometric_payoff)^2);

elapsed_time = toc; 

function result=getg0(x,kappa,theta)
    result=ones(size(x));
    result(1,:)=0*result(1,:);
    result(2,:)=result(2,:)*kappa*theta;
end

end