M = linearisedInvader(false);
ax = .02;
ay = .015;
b = 1;
c = 1;
d = 1;
mu = .24;
%%
nump = 5000;
f = zeros(nump,nump,2);
uppercutoff = 1/3.5;
dtau = ((b-mu)/c)/(nump-1) * uppercutoff;
tau = 0:dtau:(b-mu)/c *uppercutoff;

%%
X = zeros(nump,2);
Y = zeros(nump,2);
for r = 1:nump
    [eqx,eqy,stability] = equilibriumsStability(ax,ay,b,c,d,mu,tau(r),true); %This is very slow - try to do this from (x0,0) instead of globally?
    nStab = sum(stability == -1);
    X(r,1:nStab) = eqx(stability == -1);
    Y(r,1:nStab) = eqy(stability == -1);
    
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
%%
eqs = [];
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
    taupol = tau(idx);
    fisx = dfd(1:end-1).*dfd(2:end);
    eqidx = find(fisx<0);
    eqs = [eqs tau(eqidx + 1)];
    ddfd = diff(dfd);
    stab = ddfd(eqidx);
    stability = [stability stab + 1 - 2*(stab>0)];
end

%%
pip = struct;
 pip.data = {ax,ay,b,c,d,mu,X0};
 pip.dataDesc = {'ax','ay','b','c','d','mu','X0'};
 pip.tau = tau;
 pip.f = f;
 pip.eqs = eqs;
 pip.stability = stability;
 %save('PIP.mat','pip');

% pip4.data = {ax,ay,b,c,d,mu,X0};
% pip4.dataDesc = {'ax','ay','b','c','d','mu','X0'};
% pip4.tau = tau;
% pip4.X = X;
% pip4.f = f;
% pip4.eqs = eqs;
% pip4.stability = stability;
% pips.d = .75;
% pips.ays = [.015;.0125;.015;.0125]
% pips.mus = [.24;.24;.22;.22]
% save('PIP123.mat','pip1','pip2','pip3','pip4','pips');
%%
figure(1)
imshow(((f(:,:,1)>0)*2+(f(:,:,2)>0))/3,'xdata',tau,'ydata',tau)
set(gca,'visible','on','ydir','normal')
xlabel('\tau_{r}')
ylabel('\tau_{i}')
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(eqs(stability<0),eqs(stability<0),'bo');
ueq = plot(eqs(stability>=0),eqs(stability>=0),'r-o');
legend([dummy1,dummy2,seq,ueq],{'(\partial \lambda*)/(\partial \tau) < 0','(\partial \lambda*)/(\partial \tau) > 0','stable equilibrium','unstable equilibrium'})
title('Pairwise invasibility plot')