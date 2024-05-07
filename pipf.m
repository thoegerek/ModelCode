function [f,eqs,stability] = pipf(ax,ay,b,c,d,mu,nump,uppercutoff)
%pip as a function

M = linearisedInvader(false);
f = zeros(nump,nump,2);

dtau = ((b-mu)/c)/(nump-1) * uppercutoff;
tau = 0:dtau:(b-mu)/c *uppercutoff;

X = zeros(nump,2);
Y = zeros(nump,2);
for r = 1:nump
    [eqx,eqy,stability] = equilibriumsStability(ax,ay,b,c,d,mu,tau(r),true); %This is very slow - try to do this from (x0,0) instead of globally?
    nStab = sum(stability == -1);
    X(r,1:nStab) = eqx(stability == -1);
    Y(r,1:nStab) = eqy(stability == -1);

    [Y(r,1:nStab),ord] = sort(Y(r,1:nStab));
    X(r,1:nStab) = X(r,ord);
    
    for j = 1:2
        for i = 1:nump
            [~,lambda] = eig(M(ax,ay,b,c,d,mu,tau(i),tau(r),X(r,j),Y(r,j)));
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
for j = 1:size(X,2)
    dfd = zeros(nump-1,1);
    for i = 2:nump-1
        [~,lu] = eig(M(ax,ay,b,c,d,mu,tau(i+1),tau(i),X(i,j),Y(i,j)));
        [~,ld] = eig(M(ax,ay,b,c,d,mu,tau(i-1),tau(i),X(i,j),Y(i,j)));
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
