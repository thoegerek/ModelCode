function [eqx,eqy,stability] = equilibriumsStability(ax,ay,b,c,d,mu,tau,ONLY_INTERIOR)
%finds all equilibriums (of "mymodel" for given values) and their stability
%eqx = n vector of equilibriums in x (arbitrary order)
%eqy = n vector of equilibriums in y (corresponding order)
%stability n vector of stability flags (corresponding order)
%stability flags: -1 = stable, 0 = saddle, 1 = unstable

%ONLY_INTERIOR = bool: true if excluding negative equilibriums 


syms x y 

fx = b - ax*x - mu - c*tau;
fy = b - ay*y - mu - c*tau;

D = d*tau^2*x*y;

dx = x*fx + (fy + mu)*y + D*(fx-fy);
dy = -y*mu + D*(fy-fx);

ss = solve(dx==0,dy==0,{x,y});

clear x y

eqsx = double(ss.x);
eqsy = double(ss.y);


eqx = zeros(length(eqsx),1);
eqy = zeros(length(eqsx),1);
stability = zeros(length(eqsx),1);
J = jacobian(false);
for i = 1:length(eqsx)
    if isreal(eqsx(i)) && isreal(eqsy(i)) && (~ONLY_INTERIOR || (eqsx(i)>=0 && eqsy(i)>=0))
        eqx(i) = eqsx(i);
        eqy(i) = eqsy(i);
        [~,lambdas] = eig(J(ax,ay,b,c,d,mu,tau,eqsx(i),eqsy(i)));
        if det(lambdas) < 0
            stability(i) = 0;
        elseif lambdas(1) > 0
            stability(i) = 1;
        else
            stability(i) = -1;
        end
    end
end
if ONLY_INTERIOR && length(eqx) ~= 6 %For algebraic multiplicity
   eqx = zeros(6,1);
   eqy = zeros(6,1);
end