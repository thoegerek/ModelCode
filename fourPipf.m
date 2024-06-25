function [f,eqs,stability] = fourPipf(a1,a2,a3,a4,b,c,d,mu,nump,uppercutoff)
%pip as a function

M = fourLinearisedInvader(false);
f = zeros(nump,nump,2);

dtau = ((b-mu)/c)/(nump-1) * uppercutoff;
tau = 0:dtau:(b-mu)/c *uppercutoff;

N1 = zeros(nump,2);
N2 = zeros(nump,2);
N3 = zeros(nump,2);
N4 = zeros(nump,2);
for r = 1:nump
    [eq1,eq2,eq3,eq4,stability] = fourEquilibriumsStability(a1,a2,a3,a4,b,c,d,mu,tau(r),true); %This is very slow - try to do this from (x0,0) instead of globally?
    nStab = sum(stability == -1);
    N1(r,1:nStab) = eq1(stability == -1);
    N2(r,1:nStab) = eq2(stability == -1);
    N3(r,1:nStab) = eq3(stability == -1);
    N4(r,1:nStab) = eq4(stability == -1);

%     [N2(r,1:nStab),ord] = sort(N2(r,1:nStab));
%     N1(r,1:nStab) = N1(r,ord);
    
    for j = 1:2
        for i = 1:nump
            [~,lambda] = eig(M(a1,a2,a3,a4,b,c,d,mu,N1(r,j),N2(r,j),N3(r,j),N4(r,j),tau(i),tau(r)));
            f(i,r,j) = max(diag(lambda));
        end
        if nStab == 1
            f(:,r,2) = f(:,r,1);
            break;
        end
    end

    disp([num2str(r) ' / ' num2str(nump)])
end

eqs = 0;
stability = [];
for j = 1:size(N1,2)
    dfd = zeros(nump-1,1);
    for i = 2:nump-1
        [~,lu] = eig(M(a1,a2,a3,a4,b,c,d,mu,N1(r,j),N2(r,j),N3(r,j),N4(r,j),tau(i),tau(i+1)));
        [~,ld] = eig(M(a1,a2,a3,a4,b,c,d,mu,N1(r,j),N2(r,j),N3(r,j),N4(r,j),tau(i),tau(i-1)));
        dfd(i) = max(diag(ld))-max(diag(lu));
    end
    dfd(1) = [];
    [~,idx] = min(abs(dfd));
    eqs(1) = tau(idx);
    fisx = dfd(1:end-1).*dfd(2:end);
    eqidx = find(fisx<0);
    eqs = [eqs tau(eqidx + 1)];
    ddfd = diff(dfd);
    stab = ddfd(eqidx);
    stability = [stability;stab + 1 - 2*(stab>0)];
end
