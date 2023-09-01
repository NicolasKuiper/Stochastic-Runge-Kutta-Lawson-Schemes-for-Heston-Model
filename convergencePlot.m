% Function: ConvergencePlot.m
%
%
% Purpose: Convergence plot for number of steps and number 
%          of simulations for Asian and European Options.
%
% Algorithm: Nicolas Kuiper & Martin Westberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function convergencePlot(S0,r,V0,K,K_cv,type,kappa,theta,sigma,rho,Nt,Nsim,T,R)
%% Varying step-size European
step_sizes = 50:5:450;  % Number of step sizes
option_prices = zeros(1, length(step_sizes));
%% European Monte Carlo and variance reduction techniques

figure;
subplot(3,1,1);
% Standard
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
legend({'MC','','MC AV','','MC CV',''},'Location','eastoutside','FontSize',12)
title('Monte Carlo (European Call) - Step Size Convergence','FontSize',16);

%% European Midpoint FSL and variance reduction techniques 
% Standard
%figure(2);

subplot(3,1,2);
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
for k = 1:length(step_sizes)
    Nt = step_sizes(k);  
    [~, price, ~, ~] = European_Heston_Midpoint_FSL_AV(S0,r,V0,K,T,type, ...
        kappa,theta,sigma,rho,Nt,Nsim,R);
    option_prices(k)=price;
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on
% Control 
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, price, ~, ~] = European_Heston_Midpoint_FSL_CV(S0,r,V0,K,K_cv,T, ...
        type,kappa,theta,sigma,rho,Nt,Nsim,R);
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
legend({'Midpoint','','Midpoint AV','','Midpoint CV',''},'Location','eastoutside','FontSize',12)
title('Midpoint (European Call) - Step Size Convergence','FontSize',16);

%% European Euler DSL and variance reduction techniques
% Standard
%figure(3);

subplot(3,1,3);
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL(S0,r,V0,K,T, ...
        type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(step_sizes,option_prices,'filled','bo')
hold on
plot(step_sizes, option_prices, 'b:','LineWidth',2);
hold on
% Antithetic
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL_AV(S0,r,V0,K, ...
        T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on
% Control
for k = 1:length(step_sizes)
    Nt = step_sizes(k);   
    [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL_CV(S0,r,V0,K, ...
        K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
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
legend({'Euler','','Euler AV','','Euler CV',''},'Location','eastoutside','FontSize',12)
title('Euler (European Call) - Step Size Convergence','FontSize',16);
figureSize = [100, 100, 1200, 900];          % [left, bottom, width, height]
set(gcf, 'Position', figureSize);





%% Varying step-size Asian
% 
% option_prices = zeros(1, length(step_sizes));

figure;
subplot(3,1,1);
% Standard
for k = 1:length(step_sizes)
    Nt = step_sizes(k);
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC(S0, r, V0, K, type, kappa, ...
        theta, sigma, rho, Nt, Nsim, T, R);
end
scatter(step_sizes,option_prices,'filled','bo')
hold on
plot(step_sizes, option_prices, 'b:','LineWidth',2);
hold on
% Antithetic
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC_AV(S0, r, V0, K, type, ...
        kappa, theta, sigma, rho, Nt, Nsim, T, R);
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on
% Control
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC_CV(S0, r, V0, K, type, ...
        kappa, theta, sigma, rho, Nt, Nsim, T, R);
end
scatter(step_sizes,option_prices,'filled','ro')
hold on
plot(step_sizes, option_prices, 'r:','LineWidth',2);
hold on
% Plot convergence
grid on
ylabel('Option Price','FontSize',14);
xlim([step_sizes(1), step_sizes(end)]);
legend({'MC','','MC AV','','MC CV',''},'Location','eastoutside','FontSize',12)
title('Monte Carlo (Arithmetic Asian Call) - Step Size Convergence','FontSize',16);

%% European Midpoint FSL and variance reduction techniques 
% Standard
%figure(2);

subplot(3,1,2);
for k = 1:length(step_sizes)
    Nt = step_sizes(k);
    [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL(S0,r,V0,K,T,type,kappa, ...
        theta,sigma,rho,Nt,Nsim,R);
    option_prices(k)=price;
end
scatter(step_sizes,option_prices,'filled','bo')
hold on
plot(step_sizes, option_prices, 'b:','LineWidth',2);
hold on
% Antithetic
for k = 1:length(step_sizes)
    Nt = step_sizes(k);  
    [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL_AV(S0,r,V0,K,T,type,kappa, ...
        theta,sigma,rho,Nt,Nsim,R);
    option_prices(k)=price;
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on
% Control 
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL_CV(S0,r,V0,K,T,type, ...
        kappa,theta,sigma,rho,Nt,Nsim,R);
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
legend({'Midpoint','','Midpoint AV','','Midpoint CV',''},'Location','eastoutside','FontSize',12)
title('Midpoint (Arithmetic Asian Call) - Step Size Convergence','FontSize',16);

%% European Euler DSL and variance reduction techniques
% Standard
%figure(3);

subplot(3,1,3);
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL(S0,r,V0,K,T,type, ...
        kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(step_sizes,option_prices,'filled','bo')
hold on
plot(step_sizes, option_prices, 'b:','LineWidth',2);
hold on
% Antithetic
for k = 1:length(step_sizes)
    Nt = step_sizes(k);    
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL_AV(S0,r,V0,K,T ...
        ,type,kappa,theta,sigma,rho,Nt,Nsim,R);
end
scatter(step_sizes,option_prices,'filled','go')
hold on
plot(step_sizes, option_prices, 'g:','LineWidth',2);
hold on
% Control
for k = 1:length(step_sizes)
    Nt = step_sizes(k);   
    [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL_CV(S0,r,V0,K,T, ...
         type,kappa,theta,sigma,rho,Nt,Nsim,R);
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
legend({'Euler','','Euler AV','','Euler CV',''},'Location','eastoutside','FontSize',12)
title('Euler (Arithmetic Asian Call) - Step Size Convergence','FontSize',16);
figureSize = [100, 100, 1200, 900];          % [left, bottom, width, height]
set(gcf, 'Position', figureSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Varying step-size Asian
% % 
% % figure;
% % 
% % % Asian Monte Carlo
% % for k = 1:length(Nsim_values)
% %     Nsim = Nsim_values(k);
% %     
% %     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
% % end
% % plot(Nsim_values, option_prices, 'b-x');
% % hold on
% % 
% % % Asian Antithetic Monte Carlo
% % for k = 1:length(Nsim_values)
% %     Nsim = Nsim_values(k);
% %     
% %     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC_AV(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
% % end
% % plot(Nsim_values, option_prices, 'r-x');
% % hold on
% % 
% % % Asian Control Variate Monte Carlo
% % for k = 1:length(Nsim_values)
% %     Nsim = Nsim_values(k);
% %     
% %     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC_CV(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
% % end
% % plot(Nsim_values, option_prices, 'g-x');
% % hold on
% %% 
% % Asian Midpoint FSL
% 
% 
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
% end
% 
% figure(4);
% plot(Nsim_values, option_prices, 'k');
% hold on
% 
% % Asian Antithetic Midpoint FSL 
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
% end
% plot(Nsim_values, option_prices, 'b');
% hold on
% 
% % Asian Control Variate Midpoint FSL 
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, price, ~, ~, ~, ~] =  Asian_Heston_Midpoint_FSL_CV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
% end
% plot(Nsim_values, option_prices, 'c');
% xlabel('Number of Simulations','FontSize',18);
% ylabel('Option Price','FontSize',18);
% xlim([Nsim_values(1), Nsim_values(end)]);
% legend({'Midpoint AV ','Midpoint CV'},'Location','southeast','FontSize',14)
% title('Asian Call - Midpoint FSL Simuluations Convergence','FontSize',20);
% figureSize = [100, 100, 1000, 800];          % [left, bottom, width, height]
% set(gcf, 'Position', figureSize);
% grid on
% 
% %%
% 
% % Asian Monte Euler DSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
% end
% 
% figure(5)
% plot(Nsim_values, option_prices, 'k');
% hold on
% 
% % Asian Antithetic Euler DSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
% end
% plot(Nsim_values, option_prices, 'b');
% hold on
% 
% % Asian Control Variate Euler DSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] =  Asian_Heston_Euler_DSL_CV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
% end
% plot(Nsim_values, option_prices, 'c');
% xlabel('Number of Simulations','FontSize',18);
% ylabel('Option Price','FontSize',18);
% xlim([Nsim_values(1), Nsim_values(end)]);
% legend({'Euler AV','Euler CV'},'Location','southeast','FontSize',14)
% title('Arithmetic Asian Call - Simulations Convergence','FontSize',20);
% figureSize = [100, 100, 1000, 800];          % [left, bottom, width, height]
% set(gcf, 'Position', figureSize);
% grid on
% 
