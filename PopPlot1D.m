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

nump1 = 1000;
uppercutoff = 1/3;
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
    X0s = [xPops(i),yPops(i)+1e-6];
    
    disp(['i = ' num2str(i) ' / ' num2str(nump1)])
end

[~,monind] = max(TotPopb);
mon = tau(monind);
pol = max(eq(st==-1));
[~,polind] = min(abs(pol-tau));
%%
pp1d.data = {ax,ay,b,c,d,mu,X0};
pp1d.dataDesc = {'ax','ay','b','c','d','mu','X0'};
pp1d.tau = tau;
pp1d.xPopb = xPopb;
pp1d.xPops = xPops;
pp1d.yPopb = yPopb;
pp1d.yPops = yPops;
pp1d.pol = pol;
pp1d.polind = polind;
pp1d.mon = mon;
pp1d.monind = monind;
pp1d.TotPopb = TotPopb;
pp1d.TotPops = TotPops;
%save('Pop_Plot_1D.mat','pp1d');
%%
figure(1)
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
legend({'Total population','Collapsed migration','Chosen sociality','Ideal sociality'},'location','nw')
