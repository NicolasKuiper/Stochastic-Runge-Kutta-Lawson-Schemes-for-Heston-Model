clc, clear
% Heston model parameters
S0 = 80;                        % initial stock price
r = 0.05;                       % risk-free interest rate
V0 = 0.04;                      % initial volatility
global kappa
kappa = 1.0;                  % mean reversion speed (higher kappa, faster we return to theta)
global theta
theta = 0.05;                   % long-term volatility 
global sigma
sigma = 0.2;                    % volatility of volatility
rho = -0.7;                     % correlation between stock and volatility

% Option contract parameters
K = 85;                         % strike price
K_cv = 90;                      % strike price as control variate
T = 1;                          % time to maturity
type = 'call';                  % type of contract (change to 'put' for put options)

% Time discretization
Nt = 252;                    % number of time steps
h = T/Nt;                    % time step size
Nsim = 1000000;               % number of simulations for Monte Carlo
R = exp(-r*T);               % risk-free discount rate

%% Presenting Results

% Comparison Table
%comparisonTable(S0,r,V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,R)

% Heatmap for sensitivity analysis
%sensitivityAnalysis(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R)

% Number of step size convergence 
%StepSizeConvergence(S0,r,V0,K,K_cv,type,kappa,theta,sigma,rho,Nsim,T,R)

% Number of simulations convergence 
%SimulationsConvergence(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,T,R)

% Variability observation
VariabilityObservation(S0,r,V0,K,type,kappa,theta,...
    sigma,rho,Nt,Nsim,T,R)
