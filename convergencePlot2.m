function convergencePlot(S0,r,V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R)


%% Varying step-size European
step_sizes = 10:20:450;  % Number of step sizes
step_sizes_2 = linspace(252,252*length(step_sizes),length(step_sizes));
option_prices = zeros(1, length(step_sizes));

brown = [0.5, 0.25, 0];
orange = [1, 0.5, 0];
purple = [0.5, 0, 0.5];
% 
% figure;
% 
% % European Monte Carlo
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~] = European_Heston_MC(S0, r, V0, ...
%         K, type, kappa, theta, sigma, rho, Nt, Nsim, h, R);
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'b-x');
% hold on
% 
% % European Antithetic Monte Carlo
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~] = European_Heston_MC_AV(S0,r, ...
%     V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'r-x');
% hold on
% 
% % European Control Variate Monte Carlo
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~] = European_Heston_MC_CV(S0,r, ...
%     V0,K,K_cv,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'g-x');
% hold on
% 
% % European Midpoint FSL
% for k = 1:length(step_sizes_2)
%     Nt = step_sizes_2(k);
%     
%     [~, price, ~, ~] = European_Heston_Midpoint_FSL(S0, ...
%     r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'm-o');
% hold on
% 
% % European Antithetic Midpoint FSL 
% for k = 1:length(step_sizes_2)
%     Nt = step_sizes_2(k);
%     
%     [~, price, ~, ~] = European_Heston_Midpoint_FSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'c-o');
% hold on
% 
% % European Control Variate Midpoint FSL 
% for k = 1:length(step_sizes_2)
%     Nt = step_sizes_2(k);
%     
%     [~, price, ~, ~] = European_Heston_Midpoint_FSL_CV(S0,r,V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, '-o','Color',brown);
% hold on
% 
% % European Monte Euler DSL
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'k-+');
% hold on
% 
% % European Antithetic Euler DSL
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, '-+','Color',orange);
% hold on
% 
% % European Control Variate Euler DSL
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL_CV(S0,r,V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(step_sizes, option_prices, '-+','Color',purple);
% hold on
% 
% grid on
% xlabel('Number of Step Sizes','FontSize',18);
% ylabel('Option Price','FontSize',18);
% xlim([step_sizes(1), step_sizes(end)]);
% legend({'MC','MC AV','MC CV','Midpoint','Midpoint AV','Midpoint CV','Euler','Euler AV','Euler CV'},'Location','southeast','FontSize',14)
% title('European Call - Step Size Convergence','FontSize',20);
% figureSize = [100, 100, 1000, 800];          % [left, bottom, width, height]
% set(gcf, 'Position', figureSize);
% 
% %% Varying step-size Asian
% 
% option_prices = zeros(1, length(step_sizes));
% figure;
% % Asian Monte Carlo
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'b-x');
% hold on
% 
% % Asian Antithetic Monte Carlo
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC_AV(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'r-x');
% hold on
% 
% % Asian Control Variate Monte Carlo
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC_CV(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'g-x');
% hold on
% 
% % Asian Midpoint FSL
% for k = 1:length(step_sizes_2)
%     Nt = step_sizes_2(k);
%     
%     [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'm-o');
% hold on
% 
% % Asian Antithetic Midpoint FSL 
% for k = 1:length(step_sizes_2)
%     Nt = step_sizes_2(k);
%     
%     [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, 'c-o');
% hold on
% 
% % Asian Control Variate Midpoint FSL 
% for k = 1:length(step_sizes_2)
%     Nt = step_sizes_2(k);
%     
%     [~, price, ~, ~, ~, ~] =  Asian_Heston_Midpoint_FSL_CV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, '-o','Color',brown);
% hold on
% 
% % Asian Monte Euler DSL
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price);                      
% end
% plot(step_sizes, option_prices, 'k-+');
% hold on
% 
% % Asian Antithetic Euler DSL
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, '-+','Color',orange);
% hold on
% 
% % Asian Control Variate Euler DSL
% for k = 1:length(step_sizes)
%     Nt = step_sizes(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] =  Asian_Heston_Euler_DSL_CV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price);
% end
% plot(step_sizes, option_prices, '-+','Color',purple);
% hold on
% 
% grid on
% xlabel('Number of Step Sizes','FontSize',18);
% ylabel('Option Price','FontSize',18);
% xlim([step_sizes(1), step_sizes(end)]);
% legend({'MC','MC AV','MC CV','Midpoint','Midpoint AV','Midpoint CV','Euler','Euler AV','Euler CV'},'Location','southeast','FontSize',14)
% title('Arithmetic Asian Call - Step Size Convergence','FontSize',20);
% figureSize = [100, 100, 1000, 800];          % [left, bottom, width, height]
% set(gcf, 'Position', figureSize);


%% Varying number of simulations European

Nsim_values = 1000:10000:1000000;  % Number of simulations
option_prices = zeros(1, length(Nsim_values));

relErr = @(x,y) abs((x - y)) / y;
%% MC European
[~, BS_price] = BS_option_price(S0,K,sigma,r,T,type);
BS_price = repmat(BS_price,[1,length(Nsim_values)]);
% BS benchmark


figure;
% European Monte Carlo
for k = 1:length(Nsim_values)
    Nsim = Nsim_values(k);
    
    [~, option_prices(k), ~, ~] = European_Heston_MC(S0, r, V0, ...
        K, type, kappa, theta, sigma, rho, Nt, Nsim, h, R);
    option_prices(k)=relErr(option_prices,BS_price);
end
plot(Nsim_values, option_prices, 'b');
hold on

% European Antithetic Monte Carlo
for k = 1:length(Nsim_values)
    Nsim = Nsim_values(k);
    
    [~, option_prices(k), ~, ~] = European_Heston_MC_AV(S0,r, ...
    V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
    option_prices(k)=relErr(option_prices,BS_price)
end
plot(Nsim_values, option_prices, 'r');
hold on

% European Control Variate Monte Carlo
for k = 1:length(Nsim_values)
    Nsim = Nsim_values(k);
    
    [~, option_prices(k), ~, ~] = European_Heston_MC_CV(S0,r, ...
    V0,K,K_cv,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
    option_prices(k)=relErr(option_prices,BS_price)
end
plot(Nsim_values, option_prices, 'c');

grid on
xlabel('Number of Simulations','FontSize',18);
ylabel('Option Price','FontSize',18);
xlim([Nsim_values(1), Nsim_values(end)]);
legend({'Black-Scholes','MC','MC AV','MC CV'},'Location','southeast','FontSize',14)
title('European Call - Monte Carlo Simuluation Convergence','FontSize',20);
figureSize = [100, 100, 1000, 800];          % [left, bottom, width, height]
set(gcf, 'Position', figureSize);
hold off

%% Midpoint European
% figure
% plot(Nsim_values, BS_price, 'r-o','LineWidth',2);
% hold on
% % European Midpoint FSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, price, ~, ~] = European_Heston_Midpoint_FSL(S0, ...
%     r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, 'm-o');
% hold on
% 
% figure
% % European Antithetic Midpoint FSL 
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, price, ~, ~] = European_Heston_Midpoint_FSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, 'c-o');
% hold on
% 
% % European Control Variate Midpoint FSL 
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, price, ~, ~] = European_Heston_Midpoint_FSL_CV(S0,r,V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, '-o','Color',brown);
% 
% grid on
% xlabel('Number of Simulations','FontSize',18);
% ylabel('Option Price','FontSize',18);
% xlim([Nsim_values(1), Nsim_values(end)]);
% legend({'Black-Scholes','Midpoint','Midpoint AV','Midpoint CV'},'Location','southeast','FontSize',14)
% title('European Call - Midpoint FSL Simuluation Convergence','FontSize',20);
% figureSize = [100, 100, 1000, 800];          % [left, bottom, width, height]
% set(gcf, 'Position', figureSize);
% hold off
%% Euler European
% figure
% plot(Nsim_values, BS_price, 'r-o','LineWidth',2);
% hold on
% % European Monte Euler DSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, 'k-+');
% hold on
% 
% % European Antithetic Euler DSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, '-+','Color',orange);
% hold on
% 
% % European Control Variate Euler DSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~] = European_Heston_Euler_DSL_CV(S0,r,V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, '-+','Color',purple);
% 
% 
% grid on
% xlabel('Number of Simulations','FontSize',18);
% ylabel('Option Price','FontSize',18);
% xlim([Nsim_values(1), Nsim_values(end)]);
% legend({'Black-Scholes','Euler','Euler AV','Euler CV'},'Location','southeast','FontSize',14)
% title('European Call - Euler DSL Simuluation Convergence','FontSize',20);
% figureSize = [100, 100, 1000, 800];          % [left, bottom, width, height]
% set(gcf, 'Position', figureSize);
% hold off

% %% Varying step-size Asian
% 
% figure;
% 
% % Asian Monte Carlo
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, 'b-x');
% hold on
% 
% % Asian Antithetic Monte Carlo
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC_AV(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, 'r-x');
% hold on
% 
% % Asian Control Variate Monte Carlo
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_MC_CV(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, 'g-x');
% hold on
% 
% % Asian Midpoint FSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, 'm-o');
% hold on
% 
% % Asian Antithetic Midpoint FSL 
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, price, ~, ~, ~, ~] = Asian_Heston_Midpoint_FSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, 'c-o');
% hold on
% 
% % Asian Control Variate Midpoint FSL 
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, price, ~, ~, ~, ~] =  Asian_Heston_Midpoint_FSL_CV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
%     option_prices(k)=price;
% end
% plot(Nsim_values, option_prices, '-o','Color',brown);
% hold on
% 
% % Asian Monte Euler DSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, 'k-+');
% hold on
% 
% % Asian Antithetic Euler DSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] = Asian_Heston_Euler_DSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, '-+','Color',orange);
% hold on
% 
% % Asian Control Variate Euler DSL
% for k = 1:length(Nsim_values)
%     Nsim = Nsim_values(k);
%     
%     [~, option_prices(k), ~, ~, ~, ~] =  Asian_Heston_Euler_DSL_CV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
%     option_prices(k)=relErr(option_prices,BS_price)
% end
% plot(Nsim_values, option_prices, '-+','Color',purple);
% hold on
% 
% grid on
% xlabel('Number of Simulations','FontSize',18);
% ylabel('Option Price','FontSize',18);
% xlim([step_sizes(1), step_sizes(end)]);
% legend({'MC','MC AV','MC CV','Midpoint','Midpoint AV','Midpoint CV','Euler','Euler AV','Euler CV'},'Location','southeast','FontSize',14)
% title('Arithmetic Asian Call - Simulations Convergence','FontSize',20);
% figureSize = [100, 100, 1000, 800];          % [left, bottom, width, height]
% set(gcf, 'Position', figureSize);

