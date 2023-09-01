% Function: StepSizeConvergence.m
%
%
% Purpose: Convergence plot for number of step sizes for 
%          Asian and European Options.
%
% Algorithm: Nicolas Kuiper & Martin Westberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StepSizeConvergence(S0,r,V0,K,K_cv,type,kappa,theta,...
    sigma,rho,Nsim,T,R)
%% Varying step-size European
step_sizes = 50:10:800;  % Number of step sizes

%% European Monte Carlo and variance reduction techniques

% Standard
figure;
subplot(3,1,1);
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);
    [~, option_prices(k), ~, ~] = European_Heston_MC(S0, r, V0, ...
        K, type, kappa, theta, sigma, rho, Nt, Nsim, T, R);
end
scatter(step_sizes,option_prices,'filled','bo')
hold on
plot(step_sizes, option_prices, 'b:','LineWidth',2);
hold on

% Antithetic
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~] = European_Heston_MC_AV(S0,r, ...
    V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,T,R);
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on

% Control
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~] = European_Heston_MC_CV(S0,r, ...
    V0,K,K_cv,type,kappa,theta,sigma,rho,Nt,Nsim,T,R);
end
scatter(step_sizes,option_prices,'filled','ro')
hold on
plot(step_sizes, option_prices, 'r:','LineWidth',2);
hold on

% Plot convergence
grid on
ylabel('Option Price','FontSize',14);
xlim([step_sizes(1), step_sizes(end)]);
legend({'MC','','MC AV','','MC CV',''},'Location',...
    'southeast','FontSize',9)
title('Monte Carlo (European Call) - Step Size Convergence',...
    'FontSize',16);

%% European Midpoint FSL and variance reduction techniques 

% Standard
subplot(3,1,2);
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);
    [~, price, ~, ~] = European_Heston_Midpoint_FSL(S0, ...
    r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    option_prices(k)=price;
end
scatter(step_sizes,option_prices,'filled','bo')
hold on
plot(step_sizes, option_prices, 'b:','LineWidth',2);
hold on

% Antithetic
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);  
    [~, price, ~, ~] = European_Heston_Midpoint_FSL_AV(S0,r,V0,K,...
        T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    option_prices(k)=price;
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on

% Control 
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, price, ~, ~] = European_Heston_Midpoint_FSL_CV(S0,r,V0,K,...
        K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    option_prices(k)=price;
end
scatter(step_sizes,option_prices,'filled','ro')
hold on
plot(step_sizes, option_prices, 'r:','LineWidth',2);
hold on

% Plot convergence
grid on
ylabel('Option Price','FontSize',14);
xlim([step_sizes(1), step_sizes(end)]);
legend({'Midpoint','','Midpoint AV','','Midpoint CV',''},...
    'Location','southeast','FontSize',9)
title('Midpoint (European Call) - Step Size Convergence',...
    'FontSize',16);

%% European Euler DSL and variance reduction techniques

% Standard
subplot(3,1,3);
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL(S0,r...
        ,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(step_sizes,option_prices,'filled','bo')
hold on
plot(step_sizes, option_prices, 'b:','LineWidth',2);
hold on

% Antithetic
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL_AV(S0,r,...
        V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on

% Control
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);   
    [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL_CV(S0,r,...
        V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(step_sizes,option_prices,'filled','ro')
hold on
plot(step_sizes, option_prices, 'r:','LineWidth',2);
hold on

% Plot convergence
grid on
xlabel('Number of Step Sizes','FontSize',14);
ylabel('Option Price','FontSize',14);
xlim([step_sizes(1), step_sizes(end)]);
legend({'Euler','','Euler AV','','Euler CV',''},'Location',...
    'southeast','FontSize',9)
title('Euler (European Call) - Step Size Convergence',...
    'FontSize',16);
figureSize = [50, 50, 1200, 900];          
set(gcf, 'Position', figureSize);

%% Varying step-size Asian

%% European Monte Carlo and variance reduction techniques

% Standard
figure;
subplot(3,1,1);
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC(S0, r,...
        V0, K, type, kappa, theta, sigma, rho, Nt, Nsim, T, R);
end
scatter(step_sizes,option_prices,'filled','bo')
hold on
plot(step_sizes, option_prices, 'b:','LineWidth',2);
hold on

% Antithetic
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC_AV(S0, r,...
        V0, K, type, kappa, theta, sigma, rho, Nt, Nsim, T, R);
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on

% Control
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC_CV(S0, r,...
        V0, K, type, kappa, theta, sigma, rho, Nt, Nsim, T, R);
end
scatter(step_sizes,option_prices,'filled','ro')
hold on
plot(step_sizes, option_prices, 'r:','LineWidth',2);
hold on

% Plot convergence
grid on
ylabel('Option Price','FontSize',14);
xlim([step_sizes(1), step_sizes(end)]);
legend({'MC','','MC AV','','MC CV',''},'Location','southeast',...
    'FontSize',9)
title('Monte Carlo (Arithmetic Asian Call) - Step Size Convergence',...
    'FontSize',16);

%% Asian Midpoint FSL and variance reduction techniques

% Standard
subplot(3,1,2);
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);
    [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL(S0,r,V0,K,...
        T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    option_prices(k)=price;
end
scatter(step_sizes,option_prices,'filled','bo')
hold on
plot(step_sizes, option_prices, 'b:','LineWidth',2);
hold on

% Antithetic
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);  
    [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL_AV(S0,r,V0,...
        K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    option_prices(k)=price;
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on

% Control
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL_CV(S0,r,V0,...
        K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
    option_prices(k)=price;
end
scatter(step_sizes,option_prices,'filled','ro')
hold on
plot(step_sizes, option_prices, 'r:','LineWidth',2);
hold on

% Plot convergence
grid on
ylabel('Option Price','FontSize',14);
xlim([step_sizes(1), step_sizes(end)]);
legend({'Midpoint','','Midpoint AV','','Midpoint CV',''},'Location',...
    'southeast','FontSize',9)
title('Midpoint (Arithmetic Asian Call) - Step Size Convergence',...
    'FontSize',16);

%% Asian Euler DSL and variance reduction techniques

% Standard
subplot(3,1,3);
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL(S0,r,...
        V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(step_sizes,option_prices,'filled','bo')
hold on
plot(step_sizes, option_prices, 'b:','LineWidth',2);
hold on

% Antithetic
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL_AV(S0,r,...
        V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on

% Control
option_prices = zeros(1, length(step_sizes));
for k = 1:length(step_sizes)
    Nt = step_sizes(k);   
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL_CV(S0,...
        r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(step_sizes,option_prices,'filled','ro')
hold on
plot(step_sizes, option_prices, 'r:','LineWidth',2);
hold on

% Plot convergence
grid on
xlabel('Number of Step Sizes','FontSize',14);
ylabel('Option Price','FontSize',14);
xlim([step_sizes(1), step_sizes(end)]);
legend({'Euler','','Euler AV','','Euler CV',''},'Location',...
    'southeast','FontSize',9)
title('Euler (Arithmetic Asian Call) - Step Size Convergence',...
    'FontSize',16);
figureSize = [50, 50, 1200, 900];          
set(gcf, 'Position', figureSize);

end