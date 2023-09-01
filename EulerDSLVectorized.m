function [t,S,X] = EulerDSLVectorized(X0,A,g0,g,tspan,h,dW)
%Euler-Maruyame DSL scheme of strong order 0.5 and weak order 1 to approximate the solution of 
% dX=[AX+g_0(X)]dt + sum_{m=1}^Mg{m}(X)dW{m} on the time interval tspan with step
%size h, Wiener increments dW of dimension [n-1,P] where n-1 is the number
%of time steps and P the number of paths to simulate, and initial value X0. 

% Number of stochastic integrals
M = length(g);

% Create temporal grid
t = tspan(1):h:tspan(2);
n = length(t);

% Initialize solution
X = repmat(X0, [1, size(dW{1}, 2)]);
S = zeros(n-1,size(dW{1}, 2));
S(1,:) = X0(1,:);

% Do time stepping
for i = 2:n
    % Calculate exponential matrix for each time step
    Ep = expm(A * h);
    
    % Update X using the Euler-Maruyama DSL scheme
    X = X + h * g0(X);
    for m = 1:M
        X = X + g{m}(X) .* dW{m}(i - 1, :);
        X(2,:) = max(X(2,:),0);
    end
    
    % Apply the exponential matrix
    X = Ep * X;
    X(2,i) = max(X(2,i),0);
    S(i,:) = X(1,:);
end
end
