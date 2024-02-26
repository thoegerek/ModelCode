M = linearised(false);

ax = .02;
ay = .015;
b = 1;
c = 1;
d = 1;
mu = .24;
tau = .1;

% ax = 1;
% ay = 1;
% b = 1;
% c = 1;
% d = 1;
% mu = .5;
% tau = .25;


X = (0:10:100);%(0:.25:2);
Y = (0:10:100);%(0:.25:2);

nump1 = length(X);
nump2 = length(Y);

dx = zeros(nump1,nump2);
dy = zeros(nump1,nump2);

for i = 1:nump1
    for j = 1:nump2
        partial = myModel(0,[X(i);Y(j)],ax,ay,b,c,d,mu,tau);
        dx(j,i) = partial(1);
        dy(j,i) = partial(2);
    end
end

[eqx,eqy,st] =equilibriumsStability(ax,ay,b,c,d,mu,tau,true);
steqx = eqx(st==-1);
steqy = eqy(st==-1);
usteqx = eqx(st~=-1);
usteqy = eqy(st~=-1);
%%
Trange = (0:.1:30).^2;
num = 25;
rng(1)
X0 = rand(num)*max(X);
Y0 = rand(num)*max(Y);
Xt = cell(num);
C = zeros(num);
tol = 1e-4;
for i = 1:num
    for j = 1:num
        [~,Xt{i,j}] = ode45(@myModel,Trange,[X0(i,j);Y0(i,j)],[],ax,ay,b,c,d,mu,tau);
        for k = 1:2
            if (Xt{i,j}(end,1)-steqx(k))^2 + (Xt{i,j}(end,2)-steqy(k))^2 < tol
                C(i,j) = k;
            end
        end
    end
end
[~,Xtru] = ode45(@myModel,Trange,[usteqx(1)+.001;usteqy(1)],[],ax,ay,b,c,d,mu,tau);
%%
cols = {[0, 0, 0, 0.02];
    [1, 0, 0, 0.2];
    [0, 1, 0, 0.02]};
figure(1)
quiver(X,Y,dx,dy)
hold on
for i = 1:num
    for j = 1:num
        plot(Xt{i,j}(:,1),Xt{i,j}(:,2),'-','Color',cols{C(i,j)+1})
    end
end
plot(Xtru(:,1),Xtru(:,2),'-','Color',[0, 0, 0, 0.5])
psta = plot(steqx,steqy,'bo','markersize',10,'linewidth',2);
pusta = plot(usteqx,usteqy,'ro','markersize',10,'linewidth',2);
xlabel('n_{0}')
ylabel('n_{1}')
legend([psta pusta],{'Stable eq.','Unstable eq.'},'location','nw')
xlim([-2.5 100])%[-0.1 2])
ylim([-2.5 100])%[-0.1 2])

%export_fig('C:/Users/thekn/Pictures/SimpleParameters','-png','-m5')