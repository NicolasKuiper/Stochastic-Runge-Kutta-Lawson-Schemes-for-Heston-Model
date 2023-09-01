% Function: SimulationsConvergence.m
%
%
% Purpose: Convergence plot for number of simulations for 
%          Asian and European Options.
%
% Algorithm: Nicolas Kuiper & Martin Westberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SimulationsConvergence(S0,r,V0,K,type,kappa,theta,...
    sigma,rho,Nt,T,R)
%% Varying step-size European
simulations = [500 1000 5000 10000:10000:1000000];  % Number of simulations

%% European Monte Carlo
figure;
option_prices = zeros(1, length(simulations));
for k = 1:length(simulations)
    Nsim = simulations(k);
    [~, option_prices(k), ~, ~] = European_Heston_MC(S0,r,V0,K,...
        type,kappa,theta,sigma,rho,Nt,Nsim,T,R);
end
scatter(simulations,option_prices,'filled','bo')
hold on
plot(simulations, option_prices, 'b:','LineWidth',2);
hold on

%% European Midpoint FSL
option_prices = zeros(1, length(simulations));
for k = 1:length(simulations)
    Nsim = simulations(k);
    [~, option_prices(k), ~, ~] = European_Heston_Midpoint_FSL(S0,r,V0,K,T,...
        type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(simulations,option_prices,'filled','go')
hold on
plot(simulations, option_prices, 'g:','LineWidth',2);
hold on

%% European Euler DSL
option_prices = zeros(1, length(simulations));
for k = 1:length(simulations)
    Nsim = simulations(k);    
    [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL(S0,r,...
        V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(simulations,option_prices,'filled','ro')
hold on
plot(simulations, option_prices, 'r:','LineWidth',2);
hold on

%% Black-Scholes
[~, BS_price] = BS_option_price(S0,K,sigma,r,T,type);
BS_price = repmat(BS_price,[1,length(simulations)]);
plot(simulations,BS_price,'k-.')

% Plot convergence
grid on
xlabel('Number of Simulations','FontSize',14);
ylabel('Option Price','FontSize',14);
xlim([simulations(1), simulations(end)]);
legend({'MC','','Midpoint','','Euler','','Black-Scholes'},...
    'Location','eastoutside','FontSize',9)
title('European Call - Simulations Convergence','FontSize',16);
figureSize = [100, 100, 1200, 900];     
set(gcf, 'Position', figureSize);

%% Varying step-size Asian

%% Asian Monte Carlo
figure;
option_prices = zeros(1, length(simulations));
for k = 1:length(simulations)
    Nsim = simulations(k);
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC(S0, r, V0, K,...
        type, kappa, theta, sigma, rho, Nt, Nsim, T, R);
end
scatter(simulations,option_prices,'filled','bo')
hold on
plot(simulations, option_prices,'b:','LineWidth',2);
hold on

%% Asian Midpoint FSL
option_prices = zeros(1, length(simulations));
for k = 1:length(simulations)
    Nsim = simulations(k);
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL(S0,r,V0,K,T,...
        type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(simulations,option_prices,'filled','go')
hold on
plot(simulations, option_prices, 'g:','LineWidth',2);
hold on

%% Asian Euler DSL
option_prices = zeros(1, length(simulations));
for k = 1:length(simulations)
    Nsim = simulations(k);    
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL(S0,r,...
        V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(simulations,option_prices,'filled','ro')
hold on
plot(simulations, option_prices, 'r:','LineWidth',2);
hold on

% Plot convergence
grid on
xlabel('Number of Simulations','FontSize',14);
ylabel('Option Price','FontSize',14);
xlim([simulations(1), simulations(end)]);
legend({'MC','','Midpoint','','Euler CV',''},'Location',...
    'eastoutside','FontSize',9)
title('Arithmetic Asian Call - Simulations Convergence',...
    'FontSize',16);
figureSize = [100, 100, 1200, 900];          
set(gcf, 'Position', figureSize);

end