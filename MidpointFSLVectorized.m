function [t,S,X] = MidpointFSLVectorized(X0,A,g0,B,g,tspan,h,dW,g0Jac,gJac,Bexp)
%Midpoint FSL scheme of strong order 1 to approximate the solution of 
% dX=[AX+g_0(X)]dt + [BX+g_1(X)]dW on the time interval tspan with step
%size h, Wiener increments dW of dimension [n-1,P] where n-1 is the number
%of time steps and P the number of paths to simulate, and initial value X0. 
%g0Jac and g1Jac are handles to a function implementing the Jacobian matrix 
% of g0 and g1, their output should be of dimension [d,d,P] where d is the
%dimension of X0 and P is the number of different paths as described above.
%The solution of the implicit equation is solved with a single iteration of
%Newton's method (cmp. K. Debrabant and A. Kværnø, B-series analysis of
%stochastic Runge-Kutta methods that use an iterative scheme to compute 
%their internal stage values, SIAM J. Numer. Anal. 47, no. 1 (2008/09), pp. 181–203

%number of stochastic integrals
M = length(g);

for m=1:M
    if norm(A*B{m}-B{m}*A)>1e-10
        error("A and B have to commute")
    end
end

% Create temporal grid
t = tspan(1):h:tspan(2);
n = length(t);

%Initialise solution
P = size(dW{1},2);
X = repmat(X0,[1,P]);
S = zeros(n-1,P);
S(1,:) = X0(1,:);

%Prepare for vectorized solution
nw=length(X0);
zw=repmat((1:nw)',[nw,1]);
sw=kron((1:nw)',ones(nw,1));
Zws=kron((0:nw:(P-1)*nw)',ones(nw^2,1));
Sws=Zws+repmat(sw,[P,1]);
Zws=Zws+repmat(zw,[P,1]);
clear zw sw
EshpA=sparse(Zws,Sws,reshape(repmat(expm(A*h/2),1,1,P),nw^2*P,1),nw*P,nw*P,nw^2*P);
EshpAinv=inv(EshpA);
for i = 2:n-1
    %Calculate matrix exponentials
    Eshp = EshpA;
    Eshm=EshpAinv;
    for m=1:M
        [ExpBdWh,InvExpBdWh]=Bexp{m}(dW{m}(i-1,:)/2);
        Eshp = Eshp * sparse(Zws,Sws,reshape(ExpBdWh,nw^2*P,1),nw*P,nw*P,nw^2*P);
        Eshm = Eshm * sparse(Zws,Sws,reshape(InvExpBdWh,nw^2*P,1),nw*P,nw*P,nw^2*P);
    end
    Eh = Eshp*Eshp;
    
    %Update V
    Zx=CalculateZx(X);
    VmX = - Zx\f(X);
    X=X+reshape(VmX,[nw,P]);
    X = reshape(Eh*X(:),[nw,P]);
    X(2,i) = max(X(2,i),0);
    S(i,:) = X(1,:);
end

    function erg=f(Vnew)
        argument=reshape(Eshp*(X(:)+Vnew(:))/2,[nw,P]);
        temp=h*g0(argument); 
        for mind=1:M
            temp=temp+g{mind}(argument).*dW{mind}(i-1,:);
        end
        erg=Vnew(:)-X(:) -Eshm * reshape(temp,nw*P,1);
    end

    function erg=CalculateZx(XEvalPoint)
        EvalPoint=reshape(Eshp*XEvalPoint(:),[nw,P]);
        temp=h*feval(g0Jac,EvalPoint);
        for mind=1:M
            temp=temp+bsxfun(@times,feval(gJac{mind},EvalPoint),reshape(dW{mind}(i-1,:),1,1,P));
            %temp=temp+feval(gJac{mind},EvalPoint)*reshape(dW{mind}(i-1,:),1,1,P);
        end
        erg=speye(nw*P)-Eshm*sparse(Zws,Sws,reshape(temp,nw^2*P,1),nw*P,nw*P,nw^2*P)*Eshp/2;
    end

end