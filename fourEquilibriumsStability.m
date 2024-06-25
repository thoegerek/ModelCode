function [eq1,eq2,eq3,eq4,stability] = fourEquilibriumsStability(a1,a2,a3,a4,b,c,d,mu,tau,ONLY_INTERIOR)
%finds all equilibriums (of "mymodel" for given values) and their stability
%eq1 = n vector of equilibriums in n1 (or 2,3 4) (arbitrary order)
%stability n vector of stability flags (corresponding order)
%stability flags: -1 = stable, 0 = saddle, 1 = unstable

%ONLY_INTERIOR = bool: true if excluding negative equilibriums 


syms n1 n2 n3 n4

f1 = b - a1*n1 -c*tau - mu;
f2 = b - a2*n2 -c*tau - mu;
f3 = b - a3*n3 -c*tau - mu;
f4 = b - a4*n4 -c*tau - mu;

L12 = -d*tau^2*n1*n2*(f1-f2);
L13 = -d*tau^2*n1*n3*(f1-f3);
L14 = -d*tau^2*n1*n4*(f1-f4);
L23 = -d*tau^2*n2*n3*(f2-f3);
L24 = -d*tau^2*n2*n4*(f2-f4);
L34 = -d*tau^2*n3*n4*(f3-f4);

d1 = f1*n1 + (f2 + mu)*n2 + (f3 + mu)*n3 + (f4 + mu)*n4 - L12 - L13 - L14;
d2 = -mu*n2 + L12 - L23 - L24;
d3 = -mu*n3 + L13 + L23 - L34;
d4 = -mu*n4 + L14 + L24 + L34;

ss = solve(d1==0,d2==0,d3==0,d4==0,{n1,n2,n3,n4});

clear n1 n2 n3 n4

eqs1 = double(ss.n1);
eqs2 = double(ss.n2);
eqs3 = double(ss.n3);
eqs4 = double(ss.n4);


eq1 = zeros(length(eqs1),1);
eq2 = zeros(length(eqs1),1);
eq3 = zeros(length(eqs1),1);
eq4 = zeros(length(eqs1),1);
stability = zeros(length(eqs1),1);
J = fourJacobian(false);
for i = 1:length(eqs1)
    if isreal(eqs1(i)) && isreal(eqs2(i)) &&  isreal(eqs3(i)) &&  isreal(eqs4(i)) && (~ONLY_INTERIOR || (eqs1(i)>=0 && eqs2(i)>=0 && eqs3(i)>=0 && eqs4(i)>=0))
        eq1(i) = eqs1(i);
        eq2(i) = eqs2(i);
        eq3(i) = eqs3(i);
        eq4(i) = eqs4(i);
        [~,lambdas] = eig(J(a1,a2,a3,a4,b,c,d,mu,eqs1(i),eqs2(i),eqs3(i),eqs4(i),tau));
        if all(diag(lambdas>=0))
            stability(i) = 1;
        elseif ~all(diag(lambdas<0))
            stability(i) = 0;
        else
            stability(i) = -1;
        end
    end
end