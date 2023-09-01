function presentResults(S0,r,V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R)
%% Presenting Results

% European Option values
[~, MC_price, MC_std, MC_time] = European_Heston_MC(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
[~, MC_AV_price, MC_AV_std, MC_AV_time] = European_Heston_MC_AV(S0,r,V0,K,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
[~, MC_CV_price, MC_CV_std, MC_CV_time] = European_Heston_MC_CV(S0,r,V0,K,K_cv,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);

[~, Midpoint_price, Midpoint_std, Midpoint_time] = European_Heston_Midpoint_FSL(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
[~, Midpoint_AV_price, Midpoint_AV_std, Midpoint_AV_time] = European_Heston_Midpoint_FSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);
[~, Midpoint_CV_price, Midpoint_CV_std, Midpoint_CV_time] = European_Heston_Midpoint_FSL_CV(S0,r,V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,h,R);

[~, Euler_price, Euler_std, Euler_time] = European_Heston_Euler_DSL(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
[~, Euler_AV_price, Euler_AV_std, Euler_AV_time] = European_Heston_Euler_DSL_AV(S0,r,V0,K,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);
[~, Euler_CV_price, Euler_CV_std, Euler_CV_time] = European_Heston_Euler_DSL_CV(S0,r,V0,K,K_cv,T,type,kappa,theta,sigma,rho,Nt,Nsim,R);

n_sqrt = sqrt(Nsim);


% Create an Error Comparison Table
methods = {'Standard Monte Carlo', 'Antithetic Monte Carlo', 'Control Monte Carlo', 'Midpoint', ...
    'Antithetic Midpoint', 'Control Midpoint','Euler', 'Antithetic Euler', 'Control Euler'};
option_prices = [MC_price, MC_AV_price,MC_CV_price,Midpoint_price,Midpoint_AV_price,Midpoint_CV_price,Euler_price,Euler_AV_price,Euler_CV_price];
std_errors = [MC_std/n_sqrt, MC_AV_std/n_sqrt,MC_CV_std/n_sqrt,Midpoint_std/n_sqrt,Midpoint_AV_std/n_sqrt,Midpoint_CV_std/n_sqrt,Euler_std/n_sqrt,Euler_AV_std/n_sqrt,Euler_CV_std/n_sqrt];
elapsed_times = [MC_time, MC_AV_time,MC_CV_time,Midpoint_time,Midpoint_AV_time,Midpoint_CV_time,Euler_time,Euler_AV_time,Euler_CV_time];

comparison_table = table(methods', option_prices', std_errors', elapsed_times',...
    'VariableNames', {'Method', 'Option Price', 'Standard Error', 'Elapsed Time'});
disp(comparison_table);


% Create a Convergence Plot
steps = 1:Nt;
figure;
plot(steps, MC_price, 'b-', 'LineWidth', 2);
hold on;
plot(steps, MC_AV_price, 'r-', 'LineWidth', 2);
xlabel('Number of Steps');
ylabel('Option Price');
title('Convergence Plot');
legend('Standard Monte Carlo', 'Antithetic Monte Carlo');
grid on;
hold off;



% 
% % Create a Benchmarking Comparison Chart
% metrics = {'Accuracy', 'Computational Time'};
% standard_values = [accuracy_standard, computational_time_standard];
% antithetic_values = [accuracy_antithetic, computational_time_antithetic];
% 
% figure;
% bar([standard_values; antithetic_values]);
% xticklabels(metrics);
% legend('Standard Monte Carlo', 'Antithetic Monte Carlo');
% title('Benchmarking Comparison');
% ylabel('Value');
% 
% % Create a Box-and-Whisker Plot
% errors = {errors_standard, errors_antithetic};
% 
% figure;
% boxplot(errors);
% xticklabels({'Standard Monte Carlo', 'Antithetic Monte Carlo'});
% ylabel('Errors');
% title('Box-and-Whisker Plot');
% 
% % Create a Runtime Efficiency Chart
% method_names = {'Standard Monte Carlo', 'Antithetic Monte Carlo'};
% execution_times = [elapsed_time_standard, elapsed_time_antithetic];
% 
% figure;
% bar(execution_times);
% xticklabels(method_names);
% ylabel('Elapsed Time (seconds)');
% title('Runtime Efficiency Chart');

end
