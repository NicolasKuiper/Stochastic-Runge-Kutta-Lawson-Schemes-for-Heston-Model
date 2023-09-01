% Function: sensitivityAnalysis.m
%
%
% Purpose: Generate Heatmaps for sentitivity analysis 
%          
%
% Algorithm: Nicolas Kuiper & Martin Westberg
%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sensitivityAnalysis(S0,r,V0,K,T,type,kappa,...
    theta,sigma,rho,Nt,Nsim,R)
% Save variables for restart
V = V0;
S = S0;
k = K;

%% Sensitivity Analysis with variable Volatility and Starting Price
underlying_prices = 70:5:100;
volatilities = 0.04:0.04:0.3;

%% European Monte Carlo (Volatility vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(volatilities));
for i = 1:length(underlying_prices)
    for j = 1:length(volatilities)
        S0 = underlying_prices(i);
        V0 = volatilities(j);
        [~, option_prices(i,j), ~, ~] = European_Heston_MC(S0,r,V0,...
            K,type,kappa,theta,sigma,rho,Nt,Nsim,T,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(volatilities, underlying_prices, option_prices,'Colormap',...
    parula, 'ColorbarVisible', 'on');
xlabel('Volatility');
ylabel('Underlying Price');
title('European Call - Monte Carlo (Volatility vs. Underlying Price)');

%% European Midpoint (Volatility vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(volatilities));
for i = 1:length(underlying_prices)
    for j = 1:length(volatilities)
        S0 = underlying_prices(i);
        V0 = volatilities(j);
        [~, option_prices(i,j), ~, ~] = European_Heston_Midpoint_FSL(S0,...
            r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(volatilities, underlying_prices, option_prices,'Colormap',...
    parula, 'ColorbarVisible', 'on');
xlabel('Volatility');
ylabel('Underlying Price');
title('European Call - Midpoint FSL (Volatility vs. Underlying Price)');
 
%% European Euler (Volatility vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(volatilities));
for i = 1:length(underlying_prices)
    for j = 1:length(volatilities)
        S0 = underlying_prices(i);
        V0 = volatilities(j);
        [~, option_prices(i,j), ~, ~] = European_Heston_Euler_DSL(S0,...
            r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(volatilities, underlying_prices, option_prices,'Colormap',...
    parula, 'ColorbarVisible', 'on');
xlabel('Volatility');
ylabel('Underlying Price');
title('European Call - Euler DSL (Volatility vs. Underlying Price)');

%% Asian Monte Carlo (Volatility vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(volatilities));
for i = 1:length(underlying_prices)
    for j = 1:length(volatilities)
        S0 = underlying_prices(i);
        V0 = volatilities(j);
        [~, option_prices(i,j), ~, ~] = Asian_Heston_MC(S0, r, V0,...
            K, type, kappa, theta, sigma, rho, Nt, Nsim, T, R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(volatilities, underlying_prices, option_prices,'Colormap',...
    parula, 'ColorbarVisible', 'on');
xlabel('Volatility');
ylabel('Underlying Price');
title('Asian Call - Monte Carlo (Volatility vs. Underlying Price)');

%% Asian Midpoint (Volatility vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(volatilities));
for i = 1:length(underlying_prices)
    for j = 1:length(volatilities)
        S0 = underlying_prices(i);
        V0 = volatilities(j);
        [~, option_prices(i,j), ~, ~] = Asian_Heston_Midpoint_FSL(S0,...
            r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(volatilities, underlying_prices, option_prices,'Colormap',...
    parula, 'ColorbarVisible', 'on');
xlabel('Volatility');
ylabel('Underlying Price');
title('Asian Call - Midpoint FSL (Volatility vs. Underlying Price)');

%% Asian Euler (Volatility vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(volatilities));
for i = 1:length(underlying_prices)
    for j = 1:length(volatilities)
        S0 = underlying_prices(i);
        V0 = volatilities(j);
        [~, option_prices(i,j), ~, ~] = Asian_Heston_Euler_DSL(S0,...
            r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(volatilities, underlying_prices, option_prices,'Colormap',...
    parula, 'ColorbarVisible', 'on');
xlabel('Volatility');
ylabel('Underlying Price');
title('Asian Call - Euler DSL (Volatility vs. Underlying Price)');

%% Sensitivity Analysis with variable Strike and Starting Price
V0 = V; % restart volatility

underlying_prices = 70:5:100;
strike_prices = 70:5:100;

%% European Monte Carlo (Strike vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(strike_prices));
for i = 1:length(underlying_prices)
    for j = 1:length(strike_prices)
        S0 = underlying_prices(i);
        K = strike_prices(j);
        [~, option_prices(i,j), ~, ~] = European_Heston_MC(S0,r,V0,...
            K,type,kappa,theta,sigma,rho,Nt,Nsim,T,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(strike_prices, underlying_prices, option_prices,'Colormap',...
    summer, 'ColorbarVisible', 'on');
xlabel('Strike');
ylabel('Underlying Price');
title('European Call - Monte Carlo (Strike vs. Underlying Price)');

%% European Midpoint (Strike vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(strike_prices));
for i = 1:length(underlying_prices)
    for j = 1:length(strike_prices)
        S0 = underlying_prices(i);
        K = strike_prices(j);
        [~, option_prices(i,j), ~, ~] = European_Heston_Midpoint_FSL(S0,...
            r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(strike_prices, underlying_prices, option_prices,'Colormap',...
    summer, 'ColorbarVisible', 'on');
xlabel('Strike');
ylabel('Underlying Price');
title('European Call - Midpoint FSL (Strike vs. Underlying Price)');

%% European Euler (Strike vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(strike_prices));
for i = 1:length(underlying_prices)
    for j = 1:length(strike_prices)
        S0 = underlying_prices(i);
        K = strike_prices(j);
        [~, option_prices(i,j), ~, ~] = European_Heston_Euler_DSL(S0,...
            r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(strike_prices, underlying_prices, option_prices,'Colormap',...
    summer, 'ColorbarVisible', 'on');
xlabel('Strike');
ylabel('Underlying Price');
title('European Call - Euler DSL (Strike vs. Underlying Price)');

%% Asian Monte Carlo (Strike vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(strike_prices));
for i = 1:length(underlying_prices)
    for j = 1:length(strike_prices)
        S0 = underlying_prices(i);
        K = strike_prices(j);
        [~, option_prices(i,j), ~, ~] = Asian_Heston_MC(S0, r, V0, K,...
            type, kappa, theta, sigma, rho, Nt, Nsim, T, R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(strike_prices, underlying_prices, option_prices,'Colormap',...
    summer, 'ColorbarVisible', 'on');
xlabel('Strike');
ylabel('Underlying Price');
title('Asian Call - Monte Carlo (Strike vs. Underlying Price)');

%% Asian Midpoint (Strike vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(strike_prices));
for i = 1:length(underlying_prices)
    for j = 1:length(strike_prices)
        S0 = underlying_prices(i);
        K = strike_prices(j);
        [~, option_prices(i,j), ~, ~] = Asian_Heston_Midpoint_FSL(S0,...
            r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(strike_prices, underlying_prices, option_prices,'Colormap',...
    summer, 'ColorbarVisible', 'on');
xlabel('Strike');
ylabel('Underlying Price');
title('Asian Call - Midpoint FSL (Strike vs. Underlying Price)');

%% Asian Euler (Strike vs. Underlying Price)
option_prices = zeros(length(underlying_prices), length(strike_prices));
for i = 1:length(underlying_prices)
    for j = 1:length(strike_prices)
        S0 = underlying_prices(i);
        K = strike_prices(j);
        [~, option_prices(i,j), ~, ~] = Asian_Heston_Euler_DSL(S0,r,...
            V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(strike_prices, underlying_prices, option_prices,'Colormap',...
    summer, 'ColorbarVisible', 'on');
xlabel('Strike');
ylabel('Underlying Price');
title('Asian Call - Euler DSL (Strike vs. Underlying Price)');

%% Sensitivity Analysis with variable Rho and Starting Price
S0 = S; % restart price
K = k;  % restart volatility

rho_values = [-1 -0.7 -0.3 0 0.3 0.7 1];
volatilities = 0.04:0.04:0.3;

%% European Monte Carlo (Correlation vs. Volatility)
option_prices = zeros(length(volatilities), length(rho_values));
for i = 1:length(volatilities)
    for j = 1:length(rho_values)
        V0 = volatilities(i);
        rho = rho_values(j);
        [~, option_prices(i,j), ~, ~] = European_Heston_MC(S0,r,...
            V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,T,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(rho_values, volatilities, option_prices,'Colormap',...
    autumn, 'ColorbarVisible', 'on');
xlabel('Correlation');
ylabel('Volatility');
title('European Call - Monte Carlo (Correlation vs. Volatility)');

%% European Midpoint (Correlation vs. Volatility)
option_prices = zeros(length(volatilities), length(rho_values));
for i = 1:length(volatilities)
    for j = 1:length(rho_values)
        V0 = volatilities(i);
        rho = rho_values(j);
        [~, option_prices(i,j), ~, ~] = European_Heston_Midpoint_FSL(S0,...
            r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(rho_values, volatilities, option_prices,'Colormap',...
    autumn, 'ColorbarVisible', 'on');
xlabel('Correlation');
ylabel('Volatility');
title('European Call - Midpoint FSL (Correlation vs. Volatility)');

%% European Euler (Correlation vs. Volatility)
option_prices = zeros(length(volatilities), length(rho_values));
for i = 1:length(volatilities)
    for j = 1:length(rho_values)
        V0 = volatilities(i);
        rho = rho_values(j);
        [~, option_prices(i,j), ~, ~] = European_Heston_Euler_DSL(S0,...
            r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(rho_values, volatilities, option_prices,'Colormap',...
    autumn, 'ColorbarVisible', 'on');
xlabel('Correlation');
ylabel('Volatility');
title('European Call - Euler DSL (Correlation vs. Volatility)');

%% Asian Monte Carlo (Correlation vs. Volatility)
option_prices = zeros(length(volatilities), length(rho_values));
for i = 1:length(volatilities)
    for j = 1:length(rho_values)
        V0 = volatilities(i);
        rho = rho_values(j);
        [~, option_prices(i,j), ~, ~] = Asian_Heston_MC(S0,...
            r, V0, K, type, kappa, theta, sigma, rho, Nt, Nsim, T, R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(rho_values, volatilities, option_prices,'Colormap',...
    autumn, 'ColorbarVisible', 'on');
xlabel('Correlation');
ylabel('Volatility');
title('Asian Call - Monte Carlo (Correlation vs. Volatility)');

%% Asian Midpoint (Correlation vs. Volatility)
option_prices = zeros(length(volatilities), length(rho_values));
for i = 1:length(volatilities)
    for j = 1:length(rho_values)
        V0 = volatilities(i);
        rho = rho_values(j);
        [~, option_prices(i,j), ~, ~] = Asian_Heston_Midpoint_FSL(S0,r,...
            V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(rho_values, volatilities, option_prices,'Colormap',...
    autumn, 'ColorbarVisible', 'on');
xlabel('Correlation');
ylabel('Volatility');
title('Asian Call - Midpoint FSL (Correlation vs. Volatility)');

%% Asian Euler (Correlation vs. Volatility)
option_prices = zeros(length(volatilities), length(rho_values));
for i = 1:length(volatilities)
    for j = 1:length(rho_values)
        V0 = volatilities(i);
        rho = rho_values(j);
        [~, option_prices(i,j), ~, ~] = Asian_Heston_Euler_DSL(S0,...
            r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    end
end
figure('Position', [100, 100, 800, 600]);
heatmap(rho_values, volatilities, option_prices,'Colormap',...
    autumn, 'ColorbarVisible', 'on');
xlabel('Correlation');
ylabel('Volatility');
title('Asian Call - Euler DSL (Correlation vs. Volatility)');

end