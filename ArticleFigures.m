%%

load('Invasion.mat')
load('PIP.mat')
load('PIP2.mat')
load('PIP3.mat')
load('PIP4.mat')
load('Bifurcation_tau.mat')
load('PhasePop.mat')
load('PiB_d.mat')
load('Bistability_d.mat')
load('Pop_Plot_1D.mat')
load('Pop_tau.mat')
% load('AXAY.mat')
% load('PiB_dimless.mat')
%load('Pha.mat');
% load('Hys.mat');

% Color (and names) definitions
C = struct;
C.stable = [0.3 0.7 1];
C.unstable = [0.8500 0.3250 0.0980];
C.opt = [0.85 0.35 0.85];
C.totpop = [0 0 0];
C.nonmig = [0.25 0.65 0.1];
C.mig = [0.65 0.5 0];
C.fraction = [0 0.4470 0.7410];
C.bireg = [1 .4 .4];

names = cell(13,1);
%%
fig = figure;
names{get(gcf,'number')} = 'PopInv';


ymax = 2000;
ymin = 5;
Tcat = [inv.T1;inv.T2;inv.T3;inv.T4;inv.T5;inv.T6];
Xcat = [inv.X1;inv.X2;inv.X3;inv.X4;inv.X5;inv.X6];
plt = plot(Tcat,sum(Xcat,2),inv.T1,inv.X1,inv.T2,inv.X2,inv.T3,inv.X3,inv.T4,inv.X4,inv.T5,inv.X5,inv.T6,inv.X6,'linewidth',1.5);
hold on
plt(1).LineWidth = 1;
colororder([C.totpop;
    repmat([C.nonmig;C.mig],12,1)]);
% styles = flip({'^','^','^','^','^','^'...
%     's','s','s','s',...
%     'diamond','diamond','diamond','diamond',...
%     'o','o'});

styles = flip({'-','-','-','-','-','-','-',...
'--','--','--','--',...
'-.','-.','-.','-.',...
':',':',':',':',':',':',':',':',':',':'});

% spacing = 100;
% mkp1 = interp1(inv.T1,inv.X1,1:spacing:Tcat(end));
% mkp2 = interp1(inv.T2,inv.X2,1:spacing:Tcat(end));
% mkp3 = interp1(inv.T3,inv.X3,1:spacing:Tcat(end));
% mkp4 = interp1([inv.T4;inv.T5(2:end);inv.T6(2:end)],[inv.X4;inv.X5(2:end,:);inv.X6(2:end,:)],1:spacing:Tcat(end));
% 
% plot(1:spacing:Tcat(end),mkp1,'.')
% plot(1:spacing:Tcat(end),mkp2,'.')
% plot(1:spacing:Tcat(end),mkp3,'.')
% plot(1:spacing:Tcat(end),mkp4,'.')


lines = findobj(gca,'Type','line');
for i = 1:length(lines)
    lines(i).LineStyle= styles(i);
end

plot(inv.T2(find(max(inv.X2(:,3:4),[],2)>=ymin,1,'first')),ymin,'^','color',C.nonmig);
plot(inv.T3(find(max(inv.X3(:,3:4),[],2)>=ymin,1,'first')),ymin,'s','color',C.nonmig);
plot(inv.T4(find(max(inv.X4(:,3:4),[],2)>=ymin,1,'first')),ymin,'diamond','color',C.nonmig);


nothing = plot(nan,'color',[0 0 0 0]);
dummy = plot(nan,nan,'-',nan,nan,'--',nan,nan,'-.',nan,nan,':','color','k');

patch([inv.T5(1) inv.T5(1) inv.T5(end) inv.T5(end)],[0 ymax ymax 0],[.5 .5 .5],'facealpha',.3,'edgealpha',0)

xlim([0,Tcat(end)])
ylim([ymin ymax])
lgd = legend([plt(2),plt(3),plt(1),nothing,dummy(1),dummy(2),dummy(3),dummy(4)],{'Non-migrating population','Migrating population','Total population',...
    '','$\sigma = 1$','$\sigma = 0.8$','$\sigma = 0.6$','$\sigma = 0.4$'},'location','nw');
lgd.NumColumns = 2;
xlabel('Time (years)')
ylabel('Adult population number (1000s)')


set(gcf,'Position',[550 350 825 575])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

%% PIP
figure
names{get(gcf,'number')} =  'PIP';

tiledlayout(4,3)
nexttile([3 3])
imshow(((pip.f(:,:,1)<0)*2+(pip.f(:,:,2)<0))/3,'xdata',pip.tau,'ydata',pip.tau)
set(gca,'visible','on','ydir','normal')
xlabel('Resident sociality ($\sigma_{n}$)','fontsize',12)
ylabel('Invader sociality ($\sigma_{m}$)','fontsize',12)
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(0,0,'o',pip.eqs(end),pip.eqs(end),'o','color',C.stable,'markersize',10,'linewidth',2);
%ueq = plot(pip.eqs(pip.stability>=0),pip.eqs(pip.stability>=0),'-o','color',C.unstable,'markersize',10,'linewidth',2);
%greensquare = plot(.2,.1,'square','color',C.nonmig,'linewidth',5);
legend([dummy1,dummy2,seq(1)],{'Positive invasion fitness','Negative invasion fitness','Stable equilibria',...
    },'location','nw','fontsize',12)

set(0,'defaulttextInterpreter','tex') 

text(.33,.2,'$+(\div)$','color','w')

text(.33,.45,'$\div(+)$')
text(.1095,.07,'$\div$','fontsize',16)

text(.2,.6,'$\div$','fontsize',24)
text(.2,.1,'$+$','color','w','fontsize',24)
text(.8,.3,'$+$','color','w','fontsize',24)

plot([.4 .6 .6 .8 .8 1],[.4 .4 .6 .6 .8 .8],'--','color',[.5 .5 .5])
plot(1,.8,'^','color',C.nonmig,'markerfacecolor',C.nonmig);
plot(.8,.6,'s','color',C.nonmig,'markerfacecolor',C.nonmig);
plot(.6,.4,'diamond','color',C.nonmig,'markerfacecolor',C.nonmig);

%
nexttile
imshow(((pip2.f(:,:,1)>0)*2+(pip2.f(:,:,2)>0))/3,'xdata',pip2.tau,'ydata',pip2.tau)
set(gca,'visible','on','ydir','normal')
xlabel('$\tau_{r}$')
ylabel('$\tau_{i}$')
%
nexttile
imshow(((pip3.f(:,:,1)>0)*2+(pip3.f(:,:,2)>0))/3,'xdata',pip3.tau,'ydata',pip3.tau)
set(gca,'visible','on','ydir','normal')
xlabel('$\tau_{r}$')
%
nexttile
imshow(((pip4.f(:,:,1)>0)*2+(pip4.f(:,:,2)>0))/3,'xdata',pip4.tau,'ydata',pip4.tau)
set(gca,'visible','on','ydir','normal')
xlabel('$\tau_{r}$')

%set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% Bifurcation in population eq.
figure
names{get(gcf,'number')} = 'Bifurcation_pop';

%stax1b
[stax,ord] = sort(bif.sx);
tstax = bif.tsx(ord);
stax1a = stax(tstax < .3);
tstax1a = tstax(tstax < .3);

stax0 = stax(logical((tstax >= .3 ).* (tstax <= .5)));
tstax0 = tstax(logical((tstax >= .3 ).* (tstax <= .5)));
%stax1b
tstax2a = tstax0(stax0 > 550);
stax2a = stax0(stax0 > 550);
[tstax2a,ord] = sort(tstax2a);
stax2a = stax2a(ord);

%stax2a
tstax1b = tstax0(stax0 <= 550);
stax1b = stax0(stax0 <= 550);

%stax2b
stax2b = stax(tstax > .5);
tstax2b = tstax(tstax > .5);

stax1 = [stax1a;stax1b];
tstax1 = [tstax1a,tstax1b];
stax2 = [stax2a;flip(stax2b)];
tstax2 = [tstax2a,flip(tstax2b)];

%ustax
[ustax,ord] = sort(bif.ux);
tustax = bif.tux(ord);
tustax = tustax(ustax>0);
ustax = ustax(ustax>0);
tustax(ustax < min(stax1)) = [];
ustax(ustax < min(stax1)) = []; % this is a hack removing an impossible result.. I don't know when the impossible result emerged, wasn't always there

%stay1 
[stay,ord] = sort(bif.sy);

tstay = bif.tsy(ord);
tstay1 = tstay(stay==0);
stay1 = stay(stay==0);

%stay2 
tstay2 = tstay(stay~=0);
stay2 = stay(stay~=0);
[tstay2,ord] = sort(tstay2);
stay2 = stay2(ord);

%ustay
[ustay,ord] = sort(bif.uy);
tustay = bif.tuy(ord);
tustay = tustay(ustay>0);
ustay = ustay(ustay>0);
tustay = [max(tustax),tustay(end),tustay];
ustay = [0;0;ustay];

ymax = 1400;
ymin = 0; % maybe -1?


pat = patch([min(tustay) min(tustay) max(bif.tsy(bif.sy==0)) max(bif.tsy(bif.sy==0))],[ymin-1 ymax*1.1  ymax*1.1 ymin-1],C.bireg);
pat.EdgeColor = C.bireg;
pat.FaceAlpha = 0.1;
pat.EdgeAlpha = 0.2;
hold on

% saddlex1 = [min(tustax) (max(ustax)+min(stax2a))/2];
% saddley1 = [min(tustay) (max(ustay)+min(stay2))/2];
% saddlex2 = [max(tstax1)+(tau(2)-tau(1)) min(stax1)];
% saddley2 = [max(tstay1)+(tau(2)-tau(1)) 0];
% plot(saddlex1(1),saddlex1(2),'o',saddley1(1),saddley1(2),'o',saddlex2(1),saddlex2(2),'o',saddley2(1),saddley2(2),'o')

plt = plot(tustay,ustay,'--',tstay1,stay1,tstay2,stay2,tustax,ustax,'--',tstax1,stax1,tstax2,stax2,'linewidth',1.5);
colororder([
    C.nonmig;
    C.nonmig;
    C.nonmig;
    C.mig;
    C.mig;
    C.mig]);
dummy = plot(nan,nan,'k-',nan,nan,'k--');
nothing = plot(nan,'color',[0 0 0 0]);

lgd = legend([plt(5),plt(2),nothing,dummy(1),dummy(2),pat],{'Non-migrants','Migrants','','Stable pop. equilibria','Unstable pop. equilibria','Bistable region'},'location','ne');
lgd.NumColumns = 2;

xlim([bif.tau(1) bif.tau(end)])
ylim([ymin ymax])
xlabel('Sociality ($\sigma$)')
ylabel('Populatoin size')

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% Phase Pop
figure
names{get(gcf,'number')} = 'PhasePop';

red = plot(polyshape([php.Seperatrix(:,1);php.Xspan(2);php.Xspan(1)],[php.Seperatrix(:,2);php.Yspan(1)-1;php.Yspan(1)-1]),'facecolor',[1 .7 .7],'edgealpha',0);
hold on
green = plot(polyshape([php.Seperatrix(:,1);php.Xspan(2);php.Xspan(1)],[php.Seperatrix(:,2);php.Yspan(2)+1;php.Yspan(2)+1]),'facecolor',[.7 1 .7],'edgealpha',0);
plot(php.Seperatrix(:,1),php.Seperatrix(:,2),'--k','linewidth',2)

for i = 1:length(php.XT)
    for j = 1:length(php.XT)
        plot(php.XT{i,j}(:,1),php.XT{i,j}(:,2),'-','Color',[0,0,0,.01])
    end
end

sols = plot(nan,nan,'color',[.5 .5 .5]);%plot(A1(:,1),A1(:,2),'-k');
%plot(A2(:,1),A2(:,2),'-k')
%plot(A3(:,1),A3(:,2),'-k')

quiver(php.X,php.Y,php.dx,php.dy,'linewidth',1.5, 'AutoScaleFactor',0.5)

psta = plot(php.steqx,php.steqy,'bo','markersize',10,'linewidth',2);
pusta = plot(php.usteqx,php.usteqy,'ro','markersize',10,'linewidth',2);
xlabel('Non-migratory population ($n_{0}$)')
ylabel('Migratory population ($n_{1}$)')
legend([psta pusta,sols,red,green],{'Stable population equilibrium','Saddle point (Unstable eq. on manifold)','Solutions from random initial states',...
    'Converges to zero mmigrators','Converges to non-zero migrators'},'location','ne')
xlim([php.Xspan(1) php.Xspan(2)])
ylim([php.Yspan(1)-.1 php.Yspan(2)+.1])

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% Bifurcation in strategy
figure
names{get(gcf,'number')} = 'Bifurcation_strat';

reg = bwboundaries(bis.Bi);
yreg = zeros(max(reg{2}(:,1))-min(reg{2}(:,1)),2);
regtau = zeros(length(bis.Tau),1);
for i = 1:length(reg{2})
    vals = find(reg{2}(:,1)==i);
    if ~isempty(vals)
        yreg(i,1) = min(reg{2}(vals,2));
        yreg(i,2) = max(reg{2}(vals,2));
        regtau(i) = bis.Tau(i);
    end
end
regtau(regtau==0) = [];
for i = length(yreg):-1:1
    if yreg(i,:) == 0
        yreg(i,:) = [];
    end
end

[upper,ord] = unique(yreg(:,2),'last');
tauupper = regtau(ord);


ymin = 0;
ymax = .7;
axis;
hold on
ptch = patch([pib.d(pib.stability(:,1)==1),fliplr(pib.d(upper))],[pib.eq(pib.stability(:,1)==1,1);flipud(tauupper)],C.bireg,'facealpha',.2,'edgealpha',0);
plt1 = [];
plt2 = [];
for i = size(pib.eq,2):-1:1
    plt1{i} = plot(pib.d(pib.stability(:,i)==-1),pib.eq(pib.stability(:,i)==-1,i),'-','linewidth',1.5,'color',C.stable);
    plt2{i} = plot(pib.d(pib.stability(:,i)==1),pib.eq(pib.stability(:,i)==1,i),'--','linewidth',1.5,'color',C.unstable);
end
plot([0 max(pib.d)],[0 0],'-','linewidth',1.5,'color',C.stable);
plot(bis.d(upper),tauupper,'--','linewidth',1.5,'color',C.unstable)

plot(1,0,'o','color',C.stable,'markersize',5,'linewidth',2);

%plot([pips.d pips.d],[ymin ymax],':','color',[.7 .7 .7])
%text(pips.d+.005,.04,'-figure 5 (top-left)','color',[.5 .5 .5])

plot([1 1],[ymin ymax],':','color',[.5 .5 .5])
text(1+.005,.03,'-figure 3','color',[.3 .3 .3])

xlim([0 .1])
ylim([ymin ymax])
xlabel('Learning rate coefficient (c)')
ylabel('Sociality ($\sigma$)')


xzoom = [.005 .025];
yzoom = [.3 .5];
rectangle('Position', [xzoom(1), yzoom(1), xzoom(2)-xzoom(1),  yzoom(2)-yzoom(1)]);
plot([xzoom(2) .08],[.45 .45],'k-')
legend([plt1{2},plt2{1},ptch],{'Stable strategy equlibrium','Unstable strategy points','Bistable region'},'location','se')
ax1 = gca;

ax2 = axes('Position',[.6 .6 .3 .3]);
    ptch = patch([pib.d(pib.stability(:,1)==1),fliplr(pib.d(upper))],[pib.eq(pib.stability(:,1)==1,1);flipud(tauupper)],C.bireg,'facealpha',.2,'edgealpha',0);
    hold on
    for i = size(pib.eq,2):-1:1
        plot(pib.d(pib.stability(:,i)==-1),pib.eq(pib.stability(:,i)==-1,i),'-','linewidth',1.5,'color',C.stable);
        plot(pib.d(pib.stability(:,i)==1),pib.eq(pib.stability(:,i)==1,i),'--','linewidth',1.5,'color',C.unstable);
    end
    plot([0 max(pib.d)],[0 0],'-','linewidth',1.5,'color',C.stable);
    plot(bis.d(upper),tauupper,'--','linewidth',1.5,'color',C.unstable)
    plot([1 1],[ymin ymax],':','color',[.7 .7 .7])
    plot(1,pip.eqs(pip.stability<0),'o','color',C.stable,'markersize',5,'linewidth',2);
    plot([1 1],pip.eqs(pip.stability>=0),'o-','color',C.unstable,'markersize',5,'linewidth',2);
    rectangle('Position', [xzoom(1), yzoom(1), xzoom(2)-xzoom(1),  yzoom(2)-yzoom(1)]);
    xlim(xzoom)
    ylim(yzoom)

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% Phase plot of critical values in strat
% figure 
% names{get(gcf,'number')} = 'CriticalVals';
% 
% contourf(1./pha.B,log(pha.D),log(Ay(pha.ind)))
% ylabel('Learning rate coefficient nat. log. (ln(c))')
% xlabel('Natural mortality rate ($\mu$)')
% set(gca,'ydir','normal')
% cb = colorbar;
% yl = ylabel(cb,'Maximum sustained migration value of CNDD in distant habitat (d_{1})','FontSize',12,'Rotation',270);
% yl.Position(1) = yl.Position(1) + 1.5;
%     
% set(gcf,'Position',[550 350 700 500])
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','interpreter'),'interpreter','latex')
    %% Pip of 4 scenarios
% fig = figure;
%names{get(gcf,'number')} = 'PIP123';

% pip123 = {pip1,pip2;pip3,pip4};
% for i = 1:2
%     for j = 1:2
%         ind = (i-1)*2+j;
%     subplot(2,2,ind)
%     imshow(((pip123{i,j}.f(:,:,1)<0)*2+(pip123{i,j}.f(:,:,2)<0))/3,'xdata',pip123{i,j}.tau,'ydata',pip123{i,j}.tau,'InitialMagnification',10);
%     set(gca,'visible','on','ydir','normal')
%     text(.01,.205,['\mu = ' num2str(pips.mus(ind))])
%     text(.01,.185,['d_{1} = ' num2str(pips.ays(ind))])
%     if i == 2
%         xlabel('\sigma_{n}')
%     end
%     if j == 1
%         ylabel('\sigma_{m}')
%     end
%     hold on
%     dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
%     dummy2 = plot(inf,inf,'squarek');
%     %seq = plot(pip123{i}.eqs(pip123{i}.stability<0),pip123{i}.eqs(pip123{i}.stability<0),'bo','markersize',5,'linewidth',1.5);
%     %ueq = plot(pip123{i}.eqs(pip123{i}.stability>=0),pip123{i}.eqs(pip123{i}.stability>=0),'ro','markersize',5,'linewidth',1.5);
%     end
% end
% lg = legend([dummy2,dummy1],{'Positive invasion fitness','Negative invasion fitness'});
% lg.Position = [.5 .5 0 0];

% set(gcf,'Position',[550 350 700 700])
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% Population for varied tau (1D)
figure
names{get(gcf,'number')} = 'tauVaried';

colororder([C.totpop;C.fraction])
yyaxis left
plot(pp1d.tau,pp1d.TotPopb,'linewidth',1.5)
hold on
plot(pp1d.tau,pp1d.TotPops,'--','linewidth',1.5)
plot(pp1d.pol,pp1d.TotPopb(pp1d.polind),'o','markersize',5,'color',C.stable,'linewidth',2)
plot(pp1d.mon,pp1d.TotPopb(pp1d.monind),'o','markersize',5,'color',C.opt,'linewidth',2)
ylabel('Population number (1000s)')
xlabel('Sociality ($\sigma$)')
xlim([0 max(pp1d.tau)])
ylim([0 pp1d.TotPopb(pp1d.monind)*1.05])
yyaxis right
plot(pp1d.tau,pp1d.yPopb./(pp1d.xPopb+pp1d.yPopb),'linewidth',1.5)
plot(pp1d.tau,pp1d.yPops./(pp1d.xPops+pp1d.yPops),'--','linewidth',1.5)
ylim([0 1])
%ylabel('Migrating fraction of population')
legend({'Total population','Non-migratory state','Chosen sociality','Ideal sociality','Migrating fraction of population'},'location','se')

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% Population for varied ay and tau
figure
names{get(gcf,'number')} = 'd1Varied';

nump1 = length(popt.tau);
nump2 = length(popt.ax);
tol = 1e-5;

surf(popt.ax,popt.tau,popt.TotPop,popt.prop,'edgecolor','none');
collapse = popt.TotPop; collapse(popt.prop>tol) = nan;
hold on
surf(popt.ax,popt.tau,collapse,popt.prop-.1,'edgecolor','none');
view(45,55)
colormap(customcolormap([0 .5 .99 1],{'#ffb000',rgb2hex(C.mig),rgb2hex(C.nonmig),rgb2hex(C.nonmig/2)}))
c = colorbar;
ylabel(c,'Migrating fraction of population')

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

plot3(ones(dashlength,1)*dashmidpoint, popt.tau(1:dashlength), popt.TotPop(1:dashlength,dashmidind),':','color',C.stable,'linewidth',2.5);

plot3(pol0ax,pol0pol,pol0tot,'color',C.stable,'linewidth',2.5);
p1 = plot3(pol1ax,pol1pol,pol1tot,'color',C.stable,'linewidth',2.5);
p2 = plot3(popt.ax,popt.mon,popt.TotPop(sub2ind([nump1 nump2],popt.monind,1:nump2))+.5,'color',C.opt,'linewidth',2.5);

xlim([popt.ax(1) popt.ax(end)])
ylim([popt.tau(1) popt.tau(end)])
xlab = xlabel('Distant habitat CNDD ($d_{1}$)');
xlab.Rotation = -atan(320/190)*180/(2*pi);
xlab.Position(1) = .00055;
xlab.Position(2) = -.03;
ylab = ylabel('Sociality ($\sigma$)');
ylab.Rotation = atan(320/190)*180/(2*pi);
%ylab.Position(1) = -0.03;
ylab.Position(2) = .3;
zlabel('Total population number (1000s)')

legend([p1,p2],{'Chosen strategy','Ideal strategy'},'location','ne')

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% pop varying ax and ay
% figure
% names{get(gcf,'number')} = 'd0d1Varied';
% 
% nump1 = length(axay.ax);
% nump2 = length(axay.ay);
% 
% surf(axay.ax,axay.ay,log(axay.TotPop),log(axay.prop+1),'edgecolor','none');
% view(135,45)
% colormap('jet')
% c = colorbar;
% ylabel(c,'ln(d_{1}/d_{0} +1)')
% hold on
% 
% spacing1 = 3;  % better grid size
% for i = 1:spacing1:nump1
%     plot3(axay.ax, ones(nump2,1)*axay.ay(i), log(axay.TotPop(i,:)+1),'-k');
% end
% spacing2 = 10;  % better grid size
% for i = 1:spacing2:nump2
%     plot3(ones(nump1,1)*axay.ax(i), axay.ay, log(axay.TotPop(:,i)+1),'-k');
% end
% 
% xlim([axay.ax(1) axay.ax(end)])
% xlabel('d_{0}')
% ylabel('d_{1}')
% ylim([axay.ay(1) axay.ay(round(end/3))])
% zlabel('ln(n_{0}* + n_{1}* +1)')

%set(gcf,'Position',[550 350 700 500])
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% bifurcation of strat -dimless
% figure
% names{get(gcf,'number')} = 'Bif_strat_dimless';
% 
% axis;
% hold on
% for i = 1:size(pidi.stability,2)
%     plot(pidi.d(pidi.stability(:,i)==1),pidi.eq(pidi.stability(:,i)==1,i),'m.')
%     plot(pidi.d(pidi.stability(:,i)==-1),pidi.eq(pidi.stability(:,i)==-1,i),'k.')
% end
% xlim([0 max(pidi.d)])
% xlabel('c \cdot b^{3} \cdot \lambda^{-2} \cdot d_{0}^{-1}');
% ylabel('\sigma [b/\lambda]')
% 
% 
% 
% %% Hysteresis example
% figure
% names{get(gcf,'number')} = 'Hysteresis';
% 
% yyaxis left
% plot(1:length(hys.D),hys.D,'k',1:length(hys.D),hys.tau,'--b')
% yyaxis right
% plot(1:length(hys.D),hys.eq)
% xlabel('t [years]')
% ylabel('Individuals \cdot 10^{3}')
% legend('c','\sigma*','n_{0}','n_{1}')
% 
% %% One invasion
% figure
% names{get(gcf,'number')} = 'Invasion';
% 
% plot(inv.T,inv.X(:,1),inv.T,inv.X(:,2),inv.T,inv.X(:,3),':',inv.T,inv.X(:,4),':','linewidth',1.5);
% set(gca,'colororder',[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[0.8500 0.3250 0.0980]])
% hold on
% plot(inv.T,sum(inv.X,2),'k--','linewidth',1.5)
% legend('n_{1} (\sigma = 1.25)','n_{2} (\sigma = 1.25)','m_{1} (\sigma = 0.95)','m_{2} (\sigma = 0.95)','Total population')
% ylabel('Adult individuals')
% xlabel('Time [years]')

%set(gcf,'Position',[550 350 700 500])
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% exporting

path = 'C:/Users/thekn/Pictures/Article1/';

for i = 1:length(names)
    if ishandle(i)
        figure(i)
        export_fig([path,names{i}],'-png','-transparent','-m5')
    end
end
