function [eqs,stability] = evolutionaryEq(model,X0,ax,ay,b,c,d,mu,tau,tol)
%Numerical solve of polyOpt (only checks tau and outputs the best option)
%NOTE: only works with myModel right now, as "M" is staticly defined
%4th and 5th args in runToSS determines precision
%Model = @myModel
%tol = number: discard solutions where invasion fitness changes less that this number

%eqs = vector (unknown length) of equilibrium points for parameters
%stability = vector (same length) of stability flags:
%stability flags: -1 = stable, 0 = saddle, 1 = unstable


M = linearisedInvader(false);
dfd = zeros(length(tau)-1,1);

for i = 2:length(tau)-1
    [~,X,SUCCESS] = runToSS(model,1,X0,10,1e-3,{ax,ay,b,c,d,mu,tau(i)}); %Find pop. eq. by ode
    if ~SUCCESS %when MaxIt does not run to tolerance (in 10), use analytical solution instead
        Eq = zeros(6,2);
        [Eq(1:6,1),Eq(1:6,2),St] = equilibriumsStability(ax,ay,b,c,d,mu,tau(i),true);
        if length(St==-1) == 1
            X = Eq(St==-1);
        else %Repeat process (once) if needed
            [~,X,SUCCESS] = runToSS(model,1,X(end,:),20,1e-2,{ax,ay,b,c,d,mu,tau(i)});
            if ~SUCCESS
                [~,idx] = convergedEq(Eq(:,1),St,X(end,1)); %If all fails, find equilibria from 1-D convergence (not robust for 2-D)
                X = [X(idx,1) X(idx,2)];
            end
        end
    end
    [~,lu] = eig(M(ax,ay,b,c,d,mu,tau(i+1),tau(i),X(end,1),X(end,2)));
    [~,ld] = eig(M(ax,ay,b,c,d,mu,tau(i-1),tau(i),X(end,1),X(end,2)));
    dfd(i) = max(diag(ld))-max(diag(lu)); %Difference in invasion fitness above and below the diagonal
end
dfd(1) = [];
fisx = dfd(1:end-1).*dfd(2:end);
eqidx = find(fisx<-tol);
eqs = tau(eqidx + 1); %set equilibria to where dfd changes sign (more significantly that "tol")
ddfd = diff(dfd);
stability = ddfd(eqidx);
stability = 1 - 2*(stability>0); %stability found from sign change direction
