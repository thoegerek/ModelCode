M = linearisedInvader(false);
ax = .01;
ay = .001;
b = 1;
c = .125;
d = .125;
mu = .5;

% ax = .01;
% ay = .01;
% b = 1;
% c = .125;
% d = .125;
% mu = .5;

X0 = [[(b-mu)/ax,(b-mu)/ay];
    [(b-mu)/ax,1]];
O = size(X0,2);

%[~,X] = ode45(@myModel,[0 1000],X0,[],ax,ay,b,c,d,mu,0.25);
%%
nump = 200;
f = zeros(nump,O);
flow = zeros(nump);
uppercutoff = 1/3;
dtau = ((b-mu)/c)/(nump-1) * uppercutoff;
tau = 0:dtau:(b-mu)/c *uppercutoff;

X = cell(nump,O);
%%
for o = 1:O
    for r = 1:nump
        [~,X{r,o}] = runToSS(@myModel,1,X0(o,:),1e2,1e-10,{ax,ay,b,c,d,mu,tau(r)});
        for i = 1:nump
            [~,lambda] = eig(M(ax,ay,b,c,d,mu,tau(i),tau(r),X{r,o}(end,1),X{r,o}(end,2)));
            f(i,r,o) = max(diag(lambda));
        end
        disp([num2str(r+(o-1)*nump) ' / ' num2str(nump*O)])
    end
end
%%
eqs = [];
stability = [];
for o = 1:O
dfd = zeros(nump-1,1);
for i = 2:nump-1
    [~,lu] = eig(M(ax,ay,b,c,d,mu,tau(i+1),tau(i),X{i,o}(end,1),X{i,o}(end,2)));
    [~,ld] = eig(M(ax,ay,b,c,d,mu,tau(i-1),tau(i),X{i,o}(end,1),X{i,o}(end,2)));
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
% pip1.data = {ax,ay,b,c,d,mu,X0};
% pip1.dataDesc = {'ax','ay','b','c','d','mu','X0'};
% pip1.tau = tau;
% pip1.X = X;
% pip1.f = f;
% pip1.eqs = eqs;
% pip1.stability = stability;
% ds = [.06;.1;.15];
% save('PIP123.mat','pip1','pip2','pip3','ds');
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