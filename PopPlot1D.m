M = linearisedInvader(false);

% ax = 0.01;
% ay = 0.01;
% b = 1;
% c = .1;
% d = .1;
% mu = .5;

ax = .02;
ay = .0175;
b = 1;
c = 1;
d = 1;
mu = .2;

nump1 = 1000;
uppercutoff = 1/4;
dtau = ((b-mu)/c)/(nump1-1) *uppercutoff;
tau = 0:dtau:(b-mu)/c  *uppercutoff;


tauext = [tau,dtau];

%%
Xb = cell(nump1,1);
TotPopb = zeros(nump1,1);
xPopb = zeros(nump1,1);
yPopb = zeros(nump1,1);

Xs = cell(nump1,1);
TotPops = zeros(nump1,1);
xPops = zeros(nump1,1);
yPops = zeros(nump1,1);


X0b = [(b-mu)/ax,(b-mu)/ay];
X0s = [1,1e-10];
[eq,st] = evolutionaryEq(@myModel,X0b,ax,ay,b,c,d,mu,tau,0);

for i = 1:nump1
    [~,Xb{i}] = runToSS(@myModel,1,X0b,1e3,1e-7,{ax,ay,b,c,d,mu,tau(i)});
    TotPopb(i) = sum(Xb{i}(end,:));
    xPopb(i) = Xb{i}(end,1);
    yPopb(i) = Xb{i}(end,2);

    [~,Xs{i}] = runToSS(@myModel,1,X0s,1e7,1e-7,{ax,ay,b,c,d,mu,tau(i)});
    TotPops(i) = sum(Xs{i}(end,:));
    xPops(i) = Xs{i}(end,1);
    yPops(i) = Xs{i}(end,2);
    X0s = [xPops(i),yPops(i)];
    
    disp(['i = ' num2str(i) ' / ' num2str(nump1)])
end

[~,monind] = max(TotPopb);
mon = tau(monind);
pol = max(eq(st==-1));
[~,polind] = min(abs(pol-tau));
%%
% popt.data = {ay,b,c,d,mu,X0};
% popt.dataDesc = {'ay','b','c','d','mu','X0'};
% popt.tau = tau;
% popt.ax = Ay;
% popt.X = X;
% popt.pol = pol;
% popt.polind = polind;
% popt.mon = mon;
% popt.monind = monind;
% popt.TotPop = TotPop;
% popt.prop = prop;
% save('Pop_tau.mat','popt');
%%
close(1)
figure(1)
yyaxis left
plot(tau,TotPopb,'b','linewidth',1.5)
hold on
plot(tau,TotPops,'b--','linewidth',1.5)
plot(pol,TotPopb(polind),'r.','markersize',10)
plot(mon,TotPopb(monind),'m.','markersize',10)
ylabel('Population (n_{0} + n_{1})')
xlabel('Sociality (\sigma)')
xlim([0 max(tau)])
ylim([0 TotPopb(monind)*1.05])
%plot(tau,xPopb,'g',tau,xPops,'--g')
%plot(tau,yPopb,tau,yPops,'--','color',[.7,.9,0])
yyaxis right
plot(tau,yPopb./(xPopb+yPopb))
plot(tau,yPops./(xPops+yPops),'--')
ylim([0 1])
ylabel('Migrating fraction of population')
legend({'Total population','Chosen sociality','Ideal sociality','Collapsed migration'},'location','se')
title(['d_{1} = ' num2str(Ayts)])
