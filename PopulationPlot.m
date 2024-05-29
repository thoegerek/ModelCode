M = linearisedInvader(false);

% ax = .02;
% ay = .015;
% b = 1;
% c = 1;
% d = 1;
% mu = .24;

ax = .001;
ay = .0003;
b = 1;
c = .2;
d = .015;
mu = .5;

nump2 = 1000;
maxAy = .001;
minAy = .0002;
Ay = linspace(minAy,maxAy,nump2);

nump1 = 1000;
uppercutoff = 1/2.5;
dtau = ((b-mu)/c)/(nump1-1) *uppercutoff;
tau = 0:dtau:(b-mu)/c  *uppercutoff;


tauext = [tau,dtau];

%%
X = cell(nump1,nump2);
TotPop = zeros(nump1,nump2);
xPop = zeros(nump1,nump2);
yPop = zeros(nump1,nump2);

Eq = zeros(nump2,1);
St = zeros(nump2,1);

for i = 1:nump2
    X0 = [(b-mu)/ax,(b-mu)/Ay(i)];
    [eq,st] = evolutionaryEq(@myModel,X0,ax,Ay(i),b,c,d,mu,tau,0);
    if length(eq) > size(Eq,2)
        Eq = [Eq,zeros(nump2,1)];
        St = [St,zeros(nump2,1)];
    end
    Eq(i,1:length(eq)) = eq;
    St(i,1:length(st)) = st;
    for j = 1:nump1
        [~,X{i,j}] = runToSS(@myModel,1,X0,1e3,1e-7,{ax,Ay(i),b,c,d,mu,tau(j)});

        %[~,lambda] = eig(M(Ax(i),ay,b,c,d,mu,tauext(j+1),tau(j),X{i,j}(end,1),X{i,j}(end,2)));  %Finds the invasion strength just above the diagonal
        %f(j) = max(diag(lambda));

        TotPop(j,i) = sum(X{i,j}(end,:));
        xPop(j,i) = X{i,j}(end,1);
        yPop(j,i) = X{i,j}(end,2);
    end
    disp([num2str(i) ' / ' num2str(nump2)])
end

[~,monind] = max(TotPop,[],1);
mon = tau(monind);

pol = max(Eq,[],2);
polind = ones(1,nump2);
for i = 1:nump2
    polind(i) = find(tau==pol(i));
end


prop = abs(yPop)./(xPop+yPop);
prop(isnan(prop)) = 0;
%%
popt = struct;
popt.data = {ay,b,c,d,mu,X0};
popt.dataDesc = {'ay','b','c','d','mu','X0'};
popt.tau = tau;
popt.ay = Ay;
popt.pol = pol;
popt.polind = polind;
popt.mon = mon;
popt.monind = monind;
popt.TotPop = TotPop;
popt.prop = prop;
save('Pop_tau.mat','popt');
%%
figure(1)
surf(Ay,tau,TotPop,prop,'edgecolor','none');
view(45,45)
colormap('jet')
colorbar
hold on

spacing1 = 100;  % better grid size
for i = [1:spacing1:nump1, nump1]
    plot3(Ay, ones(nump2,1)*tau(i), TotPop(i,:),'-k');
end
spacing2 = 10;  % better grid size
for i = [1:spacing2:nump2, nump2]
    plot3(ones(nump1,1)*Ay(i), tau, TotPop(:,i),'-k');
end

p1 = plot3(Ay,pol,TotPop(sub2ind([nump1 nump2],polind,1:nump2)),'r','linewidth',2);
p2 = plot3(Ay,mon,TotPop(sub2ind([nump1 nump2],monind,1:nump2)),'m','linewidth',2);

xlim([Ay(1) Ay(end)])
xlabel('ay')
ylabel('\tau')
ylim([tau(1) tau(end)])
zlabel('x* + y*')

legend([p1,p2],{'Chosen strategy','Ideal strategy'},'location','ne')
%%
[~,Ayind] = min(abs(Ay - ay));
figure(2)
plot(tau,TotPop(:,Ayind),'linewidth',1.5)
hold on
plot(pol(Ayind),TotPop(polind(Ayind),Ayind),'r.','markersize',10)
plot(mon(Ayind),TotPop(monind(Ayind),Ayind),'m.','markersize',10)
legend({'Total population','Chosen sociality','Ideal sociality'},'location','se')
title(['d_{1} = ' num2str(ay)])
ylabel('Population (n_{0} + n_{1})')
xlabel('Sociality (\sigma)')
xlim([0 max(tau)])
ylim([0 max(TotPop(:,Ayind))*1.2])