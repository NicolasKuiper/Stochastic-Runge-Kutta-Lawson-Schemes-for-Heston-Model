%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File: European_Heston_Euler_DSL_AV.m
%
% Purpose: Monte Carlo simulation of the Heston model by a Euler-
%          Maruyama Drift Stochastic Lawson scheme with Antithetic Variate
%          variance reduction technique to price European Options
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
function [type, option_price, std_deviation, elapsed_time] = European_Heston_Euler_DSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R)
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

[~,~, X1] = EulerDSLVectorized(X0, A, g0, g, tspan, h, dW1);
[~,~, X2] = EulerDSLVectorized(X0, A, g0, g, tspan, h, dW2);

% Calculate the price of the European option
if strcmp(type,'call')
    payoff = max(X1(1,:)-K,0);
    antithetic_payoff = max(X2(1,:)-K,0);
elseif strcmp(type,'put')
    payoff = max(K - X1(1,:),0);
    antithetic_payoff = max(K - X2(1,:),0);
end
payoff_mean = mean(payoff);
payoff_std = std(payoff);

antithetic_payoff_mean = mean(antithetic_payoff);
antithetic_payoff_std = std(antithetic_payoff);

option_price = R*0.5*(payoff_mean + antithetic_payoff_mean);
std_deviation = 0.5*sqrt(payoff_std^2+antithetic_payoff_std^2);
elapsed_time = toc;

function result=getg0(x,kappa,theta)
    result=ones(size(x));
    result(1,:)=0*result(1,:);
    result(2,:)=result(2,:)*kappa*theta;
end

end