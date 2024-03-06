ax = .02;
ay = .015;
b = 1;
c = 1;
d = 1;
mu = .24;

M = linearised(false);

Tau = .1:.001:.2;

%% This does not work!! find out why ?

tmp = zeros(length(Tau),1);

for i = 1:length(Tau)
    [eqx,eqy,stability] = equilibriumsStability(ax,ay,b,c,d,mu,Tau(i),true);
    eq1 = (eqy==0).*(eqx>0); 
    eqx1 = eqx(logical(eq1));

    lam = M(ax,ay,b,c,d,mu,.2,eqx1,0);
    [vec,val] = eig(lam);

    if -val(1,1) > val(2,2)*vec(1,2)
        tmp(i) = -1;
    else
        tmp(i) = 1;
    end
    disp(['i = ' num2str(i) ' / ' num2str(length(Tau))])
end