%%
%load('Invasion.mat')
load('PIP.mat')
 load('Bifurcation_tau.mat')
% load('PIP123.mat')
 load('PiB_d.mat')
% load('Pop_tau.mat')
% load('AXAY.mat')
% load('PiB_dimless.mat')
% load('Pha.mat');
% load('Hys.mat');

%% Pairwise invasability plot
figure(1)
imshow(((pip.f(:,:,1)<0)*2+(pip.f(:,:,2)<0))/3,'xdata',pip.tau,'ydata',pip.tau)
set(gca,'visible','on','ydir','normal')
xlabel('Sociality of prevailing population (\sigma_{n})','fontsize',12)
ylabel('Sociality of mutant population (\sigma_{m})','fontsize',12)
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(0,0,'o',pip.eqs(pip.stability<0),pip.eqs(pip.stability<0),'o','color',[0.1 0.4 0.95],'markersize',10,'linewidth',2);
ueq = plot(pip.eqs(pip.stability>=0),pip.eqs(pip.stability>=0),'-o','color',[0.8500 0.3250 0.0980],'markersize',10,'linewidth',2);
greensquare = plot(.2,.1,'square','color',[0.4660 0.6740 0.1880],'linewidth',5);
legend([dummy1,dummy2,seq(1),ueq,greensquare],{'Positive invasion fitness','Negative invasion fitness','Stable equilibria',...
    'Unstable equilibria','Case from figure 1'},'location','nw','fontsize',12)

text(.0971,.13,'+(\div)','color','w')
text(.1045,.125,'+','color','w')

text(.0972,.08,'\div(+)')
text(.1095,.07,'\div','fontsize',16)

text(.04,.11,'\div','fontsize',24)
text(.07,.03,'+','color','w','fontsize',24)
text(.17,.06,'+','color','w','fontsize',24)



%% Bifurcation in population eq.
figure(2)

%stax1 = 
%stax2 =
[ustax,ord] = sort(bif.ux);
tustax = bif.tux(ord);
tustax = tustax(ustax>0);
ustax = ustax(ustax>0);
%stay1 = 
%stay2 =
[ustay,ord] = sort(bif.uy);
tustay = bif.tuy(ord);
tustay = tustay(ustay>0);
ustay = ustay(ustay>0);
tustay = [max(tustax),tustay(end),tustay];
ustay = [0;0;ustay];

ymax = 45;

plot(tustax,ustax,'--',bif.tsx,bif.sx,'.',tustay,ustay,'--',bif.tsy,bif.sy,'.','linewidth',1.5);
colororder([
    0.4660 0.6740 0.1880;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250
    0.9290 0.6940 0.1250]);
hold on
xdum = plot(nan,nan,'-','color',[0.4660 0.6740 0.1880]);
ydum = plot(nan,nan,'-','color',[0.9290 0.6940 0.1250]);
dummy = plot(nan,nan,'k-',nan,nan,'k--',nan,nan,':');

plot([min(tustay) min(tustay)],[0 ymax],':',[max(bif.tsy(bif.sy==0)) max(bif.tsy(bif.sy==0))],[0 ymax],':','color',[0.8500 0.3250 0.0980],'linewidth',.5)
plot(pip.eqs(pip.stability<0),[0 ymax],':','color',[0.1 0.4 0.95],'linewidth',.5)

legend([xdum(1),ydum(1),dummy(1),dummy(2),dummy(3)],{'Non-migrants','Migrants','Stable pop. equilibrium','Unstable pop. equilibrium','Sociality equilibria'},'location','sw')
xlim([bif.tau(1) bif.tau(end)])
xlabel('Sociality (\sigma)')
ylabel('Populatoin size')

set(gcf,'Position',[550 350 700 500])
%% Bifurcation in strategy
figure(3)
ymin = -.01;
ymax = .15;
axis;
hold on
for i = size(pib.eq,2):-1:1
    plot(pib.d(pib.stability(:,i)==-1),pib.eq(pib.stability(:,i)==-1,i),'-','linewidth',1.5,'color',[0.1 0.4 0.95]);
    plot(pib.d(pib.stability(:,i)==1),pib.eq(pib.stability(:,i)==1,i),'--','linewidth',1.5,'color',[0.8500 0.3250 0.0980]);
end
plot([0 max(pib.d)],[0 0],'-','linewidth',1.5,'color',[0.1 0.4 0.95]);

plot([pips.d pips.d],[ymin ymax],':','color',[.7 .7 .7])
text(pips.d+.005,.04,'-figure 5 (top-left)','color',[.5 .5 .5])

plot([1 1],[ymin ymax],':','color',[.7 .7 .7])
text(1+.005,.03,'-figure 2','color',[.5 .5 .5])

xlim([0 5])
ylim([ymin ymax])
xlabel('Learning rate coefficient (c)')
ylabel('Sociality (\sigma)')
legend('Stable strategy equlibrium','Unstable strategy equilibrium')

set(gcf,'Position',[550 350 700 500])
%% Pip of 4 scenarios
fig = figure(4);
pip123 = {pip1,pip2;pip3,pip4};
for i = 1:2
    for j = 1:2
        ind = (i-1)*2+j;
    subplot(2,2,ind)
    imshow(((pip123{i,j}.f(:,:,1)<0)*2+(pip123{i,j}.f(:,:,2)<0))/3,'xdata',pip123{i,j}.tau,'ydata',pip123{i,j}.tau,'InitialMagnification',10);
    set(gca,'visible','on','ydir','normal')
    text(.01,.205,['\mu = ' num2str(pips.mus(ind))])
    text(.01,.185,['d_{1} = ' num2str(pips.ays(ind))])
    if i == 2
        xlabel('\sigma_{n}')
    end
    if j == 1
        ylabel('\sigma_{m}')
    end
    hold on
    dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
    dummy2 = plot(inf,inf,'squarek');
    %seq = plot(pip123{i}.eqs(pip123{i}.stability<0),pip123{i}.eqs(pip123{i}.stability<0),'bo','markersize',5,'linewidth',1.5);
    %ueq = plot(pip123{i}.eqs(pip123{i}.stability>=0),pip123{i}.eqs(pip123{i}.stability>=0),'ro','markersize',5,'linewidth',1.5);
    end
end
set(gcf,'Position',[550 350 700 700])
lg = legend([dummy2,dummy1],{'Positive invasion fitness','Negative invasion fitness'});
lg.Position = [.5 .5 0 0];
%% Population for varied ay and tau
nump1 = length(popt.tau);
nump2 = length(popt.ax);
figure(5)
surf(popt.ax,popt.tau,popt.TotPop,popt.prop,'edgecolor','none');
view(45,45)
colormap('jet')
c = colorbar;
ylabel(c,'{n1*}/{n0*}')
hold on

spacing1 = 100;  % better grid size
for i = [1:spacing1:nump1,nump1]
    plot3(popt.ax, ones(nump2,1)*popt.tau(i), popt.TotPop(i,:),'-k');
end
spacing2 = 10;  % better grid size
for i = [1:spacing2:nump2,nump2]
    plot3(ones(nump1,1)*popt.ax(i), popt.tau, popt.TotPop(:,i),'-k');
end

Totpol = popt.TotPop(sub2ind([nump1 nump2],popt.polind,1:nump2));

pol0pol = popt.pol(popt.pol == 0);
pol0ax = popt.ax(popt.pol == 0);
pol0tot = Totpol(popt.pol == 0);

pol1pol = popt.pol(popt.pol ~= 0);
pol1polind = popt.polind(popt.pol ~= 0);
pol1ax = popt.ax(popt.pol  ~= 0);
pol1tot = Totpol(popt.pol ~= 0);

dashmidpoint = (pol1ax(end)+pol0ax(1))/2;
dashmidind = length(pol1ax);
dashlength = pol1polind(end);

plot3(ones(dashlength,1)*dashmidpoint, popt.tau(1:dashlength), popt.TotPop(1:dashlength,dashmidind),':r','linewidth',2);

plot3(pol0ax,pol0pol,pol0tot,'r','linewidth',2);
p1 = plot3(pol1ax,pol1pol,pol1tot,'r','linewidth',2);
p2 = plot3(popt.ax,popt.mon,popt.TotPop(sub2ind([nump1 nump2],popt.monind,1:nump2)),'m','linewidth',2);

xlim([popt.ax(1) popt.ax(end)])
xlabel('d_{1}')
ylabel('\sigma')
ylim([popt.tau(1) popt.tau(end)])
zlabel('n_{0}*+ n_{1}*')

legend([p1,p2],{'Chosen strategy','Ideal strategy'},'location','ne')
%% pop varying ax and ay
nump1 = length(axay.ax);
nump2 = length(axay.ay);
figure(6)
surf(axay.ax,axay.ay,log(axay.TotPop),log(axay.prop+1),'edgecolor','none');
view(135,45)
colormap('jet')
c = colorbar;
ylabel(c,'ln(d_{1}/d_{0} +1)')
hold on

spacing1 = 3;  % better grid size
for i = 1:spacing1:nump1
    plot3(axay.ax, ones(nump2,1)*axay.ay(i), log(axay.TotPop(i,:)+1),'-k');
end
spacing2 = 10;  % better grid size
for i = 1:spacing2:nump2
    plot3(ones(nump1,1)*axay.ax(i), axay.ay, log(axay.TotPop(:,i)+1),'-k');
end

xlim([axay.ax(1) axay.ax(end)])
xlabel('d_{0}')
ylabel('d_{1}')
ylim([axay.ay(1) axay.ay(round(end/3))])
zlabel('ln(n_{0}* + n_{1}* +1)')


%% bifurcation of strat -dimless
figure(7)
axis;
hold on
for i = 1:size(pidi.stability,2)
    plot(pidi.d(pidi.stability(:,i)==1),pidi.eq(pidi.stability(:,i)==1,i),'m.')
    plot(pidi.d(pidi.stability(:,i)==-1),pidi.eq(pidi.stability(:,i)==-1,i),'k.')
end
xlim([0 max(pidi.d)])
xlabel('c \cdot b^{3} \cdot \lambda^{-2} \cdot d_{0}^{-1}');
ylabel('\sigma [b/\lambda]')


%% Phase plot of critical values in strat
figure(8)
logbif = log(pha.bif);
logbif(logbif==-inf) = max(logbif,[],'all')+1;
imagesc(logbif,'xdata',log(pha.Ay),'ydata',pha.Mu)
xlabel('ln(d_{1})')
ylabel('\mu')
set(gca,'ydir','normal')
cb = colorbar;
yl = ylabel(cb,'lower bound for critical value of c \cdot b^{3} \cdot \lambda^{-2} \cdot d_{0}^{-1}','FontSize',12,'Rotation',270);
yl.Position(1) = yl.Position(1) + 1.5;
%% Hysteresis example
figure(9)
yyaxis left
plot(1:length(hys.D),hys.D,'k',1:length(hys.D),hys.tau,'--b')
yyaxis right
plot(1:length(hys.D),hys.eq)
xlabel('t [years]')
ylabel('Individuals \cdot 10^{3}')
legend('c','\sigma*','n_{0}','n_{1}')

%% One invasion
figure(10)
plot(inv.T,inv.X(:,1),inv.T,inv.X(:,2),inv.T,inv.X(:,3),':',inv.T,inv.X(:,4),':','linewidth',1.5);
set(gca,'colororder',[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[0.8500 0.3250 0.0980]])
hold on
plot(inv.T,sum(inv.X,2),'k--','linewidth',1.5)
legend('n_{1} (\sigma = 1.25)','n_{2} (\sigma = 1.25)','m_{1} (\sigma = 0.95)','m_{2} (\sigma = 0.95)','Total population')
ylabel('Adult individuals')
xlabel('Time [years]')

%% exporting

path = 'C:/Users/thekn/Pictures/Article1/';
names = {'PIP';
    'Bifurcation_pop';
    'Bifurcation_strat';
    'PIP123';
    'd1Varied';
    'd0d1Varied';
    'Bif_strat_dimless';
    'CriticalVals';
    'Hysteresis';
    'Invasion'};

for i = 1:length(names)
    if ishandle(i)
        figure(i)
        export_fig([path,names{i}],'-png','-transparent','-m5')
    end
end
