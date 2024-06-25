%%
close all
%%
load('Invasion.mat')
load('PIP.mat')
load('PIP2.mat')
load('PIP3.mat')
load('PIP4.mat')
load('Bifurcation_tau.mat')
load('PhasePop.mat')
load('PiB_d.mat')
load('PiB_ay.mat')
load('PiB_mu.mat')
load('PopFromBif.mat')
load('PopFromBif_ay.mat')
load('PopFromBif_mu.mat')
load('Bistability_d.mat') 
load('Pop_Plot_1D.mat')
load('Pop_tau.mat')
% load('AXAY.mat')
% load('PiB_dimless.mat')
load('Pha.mat');
% load('Hys.mat');

% Color (and names) definitions
C = struct;
C.stable = [0.3 0.7 1];
C.unstable = [0.8500 0.3250 0.0980];
C.opt = [0.8500 0.3250 0.0980];%[.8 .2 .2];%[0.85 0.2 0.85];
C.totpop = [0 0 0];
C.nonmig = [0.6235    0.4000    0.7176];
C.mig = [0.25 0.65 0.1];
C.fraction = [0 0.4470 0.7410];
C.bireg = [.8 .2 .2];

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

styles = flip({'-','-.','-.','-.','-.','-.','-.',...
'-','-','-','-',...
'-.','-.','-.','-.',...
'-','-','-','-','-','-','-','-','-','-'});

lines = findobj(gca,'Type','line');
for i = 1:length(lines)
    lines(i).LineStyle= styles(i);
end

plot(225,ymin,'^','color',C.bireg,'markerfacecolor',C.bireg,'markersize',8);%inv.T2(find(max(inv.X2(:,3:4),[],2)>=ymin,1,'first'))
plot(800,ymin,'s','color',C.bireg,'markerfacecolor',C.bireg,'markersize',10);%inv.T3(find(max(inv.X3(:,3:4),[],2)>=ymin,1,'first'))
plot(1450,ymin,'diamond','color',C.bireg,'markerfacecolor',C.bireg,'markersize',9);%inv.T4(find(max(inv.X4(:,3:4),[],2)>=ymin,1,'first'))


text(25,max(inv.X1,[],'all')+40,'$\sigma = 1$')
text(540,max(inv.X2,[],'all')+40,'$\sigma = 0.8$')
text(1200,max(inv.X4,[],'all')+40,'$\sigma = 0.6$')
text(2200,max(inv.X5,[],'all')+40,'$\sigma = 0.4$')
text(2188,1850,'$\mu = 5.7 -$')

patch([inv.T5(1) inv.T5(1) inv.T5(end) inv.T5(end)],[0 ymax ymax 0],[.5 .5 .5],'facealpha',.3,'edgealpha',0)

xlim([0,Tcat(end)])
ylim([ymin ymax])
legend([plt(8),plt(9),plt(1)],{'Non-migrating population','Migrating population','Total population',...
    },'location','nw');
xlabel('Time')
ylabel('Population number')


set(gcf,'Position',[550 350 825 575])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

%% PIP
figure
names{get(gcf,'number')} =  'PIP';

tiledlayout(4,3)
nexttile([3 3])
imshow(((pip.f(:,:,1)<0)+(pip.f(:,:,2)<0)*2)/3,'xdata',pip.tau,'ydata',pip.tau)
set(gca,'visible','on','ydir','normal','box','off')
ylabel('Invader sociality ($\sigma_{m}$)\hspace{100pt}','fontsize',12)
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(0,0,'o',pip.eqs(end),pip.eqs(end),'o','color',C.stable,'markersize',10,'linewidth',2);
%ueq = plot(pip.eqs(pip.stability>=0),pip.eqs(pip.stability>=0),'-o','color',C.unstable,'markersize',10,'linewidth',2);


set(0,'defaulttextInterpreter','tex') 

text(.33,.2,'$+(-)$','color','w')

text(.33,.5,'$-(+)$')
text(.1095,.07,'$-$','fontsize',16)

text(.185,.6,'$-$','fontsize',24)
text(.185,.1,'$+$','color','w','fontsize',24)
text(.8,.3,'$+$','color','w','fontsize',24)

plot([.4 .6 .6 .8 .8 1],[.4 .4 .6 .6 .8 .8],'--','color',[.5 .5 .5],'linewidth',1.5);
plot(.5,.4,'<','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])
plot(.6,.5,'v','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])
plot(.7,.6,'<','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])
plot(.8,.7,'v','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])
plot(.9,.8,'<','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])

example = plot(nan,nan,'<--','color',[.5 .5 .5],'linewidth',1.5,'markerfacecolor',[.5 .5 .5]);

plot(1,.8,'^','color',C.bireg,'markerfacecolor',C.bireg,'markersize',8);
plot(.8,.6,'s','color',C.bireg,'markerfacecolor',C.bireg,'markersize',10);
plot(.6,.4,'diamond','color',C.bireg,'markerfacecolor',C.bireg,'markersize',9);

legend([seq(1),example],{'Evolutionary Stable Strategy','Path of adaptation (fig. 1)'...
    },'location','ne','fontsize',12);


nexttile
imshow(((pip4.f(:,:,1)<0)+(pip4.f(:,:,2)<0)*2)/3,'xdata',pip4.tau,'ydata',pip4.tau)
hold on
set(gca,'visible','on','ydir','normal','box','off')
plot(0,0,'o','color',C.stable,'markersize',7,'linewidth',1.5)
text(.025,.9,'$d_1\uparrow$')
%
nexttile
imshow(((pip2.f(:,:,1)<0)+(pip2.f(:,:,2)<0)*2)/3,'xdata',pip2.tau,'ydata',pip2.tau)
hold on
set(gca,'visible','on','ydir','normal','box','off')
plot(0,0,'o',max(pip2.eqs),max(pip2.eqs),'o','color',C.stable,'markersize',7,'linewidth',1.5)
xlabel('Resident sociality ($\sigma_{n}$)','fontsize',12)
text(.025,.9,'$d_1\uparrow$')
text(.025,.75,'$c\uparrow$')
%
nexttile
imshow(((pip3.f(:,:,1)<0)+(pip3.f(:,:,2)<0)*2)/3,'xdata',pip3.tau,'ydata',pip3.tau)
hold on
set(gca,'visible','on','ydir','normal','box','off')
plot(0,0,'o',max(pip3.eqs),max(pip3.eqs),'o','color',C.stable,'markersize',7,'linewidth',1.5)

text(.025,.9,'$c\uparrow$')
text(.025,.75,'$\lambda\downarrow$')


set(gcf,'Position',[680          91         713        1007])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
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
pat.FaceAlpha = 0.2;
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

lgd = legend([plt(5),plt(2),dummy(1),dummy(2)],{'Non-migrants','Migrants','Stable eq.','Saddle points'},'location','nw');
lgd.NumColumns = 2;

text(.02,1100,'\textbf{I}')
text(.365,1100,'\textbf{II}')
text(.935,1100,'\textbf{III}')

xlim([bif.tau(1) bif.tau(end)])
ylim([ymin ymax])
xlabel('Sociality ($\sigma$)')
ylabel('Population number')

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% Phase Pop
figure
names{get(gcf,'number')} = 'PhasePop';

plot(polyshape([php.Seperatrix(:,1);php.Xspan(2);php.Xspan(1)],[php.Seperatrix(:,2);php.Yspan(1)-1;php.Yspan(1)-1]),'facecolor',C.bireg,'edgealpha',0);
hold on
plot(polyshape([php.Seperatrix(:,1);php.Xspan(2);php.Xspan(1)],[php.Seperatrix(:,2);php.Yspan(2)+1;php.Yspan(2)+1]),'facecolor',[.8 1 .8],'edgealpha',0);
red = plot(polyshape(nan,nan),'facecolor',C.bireg,'edgecolor','k');
green = plot(polyshape(nan,nan),'facecolor',[.8 1 .8],'edgecolor','k');
plot(php.Seperatrix(:,1),php.Seperatrix(:,2),'--k','linewidth',2)

for i = 1:length(php.XT)
    for j = 1:length(php.XT)
        plot(php.XT{i,j}(:,1),php.XT{i,j}(:,2),'-','Color',[0,0,0,.01])
    end
end

sols = plot(nan,nan,'color',[.5 .5 .5]);
manif = plot(nan,nan,'k','linewidth',1.5);

norms = sqrt(php.dx.^2+php.dy.^2);
quiver(php.X,php.Y,php.dx./sqrt(norms),php.dy./sqrt(norms),'linewidth',3, 'AutoScaleFactor',1,'color',[.9 .9 .9])
quiver(php.X,php.Y,php.dx./sqrt(norms),php.dy./sqrt(norms),'linewidth',1, 'AutoScaleFactor',1,'color',[.4 .4 .4])

psta = plot(php.steqx,php.steqy,'o','color',C.stable,'markersize',10,'linewidth',2);
pusta = plot(php.usteqx(1),php.usteqy(1),'o','color',C.unstable,'markersize',10,'linewidth',2);
psad = plot(php.usteqx(2),php.usteqy(2),'o','color',[.4 .4 .4],'markersize',10,'linewidth',2);
xlabel('Non-migrant population ($n_{0}$)')
ylabel('Migrant population ($n_{1}$)')
legend([psta,psad,pusta,sols,manif,green,red],{'Stable eq.','Saddle point','Unstable eq.','Random trajectories','Attracting manifold'...
    'Basin of upper eq.','Basin of lower eq.'},'location','ne')
xlim([php.Xspan(1) php.Xspan(2)])
ylim([php.Yspan(1)-.1 php.Yspan(2)+.1])

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
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
patch([min(pib.d(pib.stability(:,1)==1)),pib.d(pib.stability(:,1)==1),fliplr(pib.d(upper))],[1;pib.eq(pib.stability(:,1)==1,1);flipud(tauupper)],C.bireg,'facealpha',.2,'edgealpha',0);
ptch = patch(nan,nan,C.bireg,'edgecolor',C.unstable,'facealpha',.2);
plt1 = [];
plt2 = [];
plot([0 max(pib.d)],[0 0],'-','linewidth',1.5,'color',C.stable);
plot(bis.d(upper),tauupper,'-','linewidth',1,'color',C.unstable)
plot([min(pib.d(pib.stability(:,1)==1)),min(pib.d(pib.stability(:,1)==1))],[max(pib.eq(pib.stability(:,i)==1,i)) tauupper(1)],'-','linewidth',1,'color',C.unstable)
for i = size(pib.eq,2):-1:1
    plt1{i} = plot(pib.d(pib.stability(:,i)==-1),pib.eq(pib.stability(:,i)==-1,i),'-','linewidth',1.5,'color',C.stable);
    plt2{i} = plot(pib.d(pib.stability(:,i)==1),pib.eq(pib.stability(:,i)==1,i),'-','linewidth',1,'color',C.unstable);
end

xlim([0 .1])
ylim([ymin ymax])
xlabel('Learning rate coefficient (c)')

ylabel('Sociality ($\sigma$)')

% qx = linspace(.005,.095,12);
% qy = linspace(.05,.75,11);
% [QX,QY] = meshgrid(qx,qy);
% 
% updown = ones(length(qy),length(qx));
% both = updown;
% updown(QX<pib.d(find(pib.eq(:,1)~=inf,1,'first'))) = deal(-1);
% qyequp = interp1(pib.d,pib.eq(:,2),qx);
% qyeqdown = interp1(pib.d,pib.eq(:,1),qx);
% qybi = interp1(bis.d(upper),tauupper,qx);
% updown(QY>qyequp) = deal(-1);
% updown(QY<qyeqdown) = deal(-1);
% updown(logical((QY>qyeqdown).*(QY<qybi))) = deal(0);
% both(logical(1-(QY>qyeqdown).*(QY<qybi))) = deal(0);
% 
% updown(1:2:end) = 0;
% both(1:2:end) = 0;

%pltq = quiver(QX,QY,zeros(length(qy),length(qx)),updown,'linewidth',1.5, 'AutoScaleFactor',.25,'color',C.quiver);
%quiver(QX,QY,zeros(length(qy),length(qx)),both,'linewidth',1.5, 'AutoScaleFactor',.25,'color',C.quiver);
%quiver(QX,QY,zeros(length(qy),length(qx)),-both,'linewidth',1.5, 'AutoScaleFactor',.25,'color',C.quiver);
grid on
set(gca,'xGrid','off','GridLineStyle','--','gridalpha',.5);

legend([plt1{2},ptch],{'ESS','Bistable region'},'location','sw')

Ls = pib.eq(pib.stability(:,2)==-1,2).^2.*pfb.D'.*pfb.X(:,1).*pfb.X(:,2).*(pfb.X(:,1)*pfb.data{1}-pfb.X(:,2)*pfb.data{2});
Lul = pib.eq(pib.stability(:,1)==1,1).^2.*pfb.D'.*pfb.X(:,1).*pfb.X(:,2).*(pfb.X(:,1)*pfb.data{1}-pfb.X(:,2)*pfb.data{2});
Luu = [nan;nan;nan;nan;tauupper].^2.*pfb.D'.*pfb.X(:,1).*pfb.X(:,2).*(pfb.X(:,1)*pfb.data{1}-pfb.X(:,2)*pfb.data{2});
Lgrid = zeros(length(Ls),6);
for i = 1:10
    Lgrid(:,i) = (i*.1)^2*pfb.D'.*pfb.X(:,1).*pfb.X(:,2).*(pfb.X(:,1)*pfb.data{1}-pfb.X(:,2)*pfb.data{2}); %works without the ^2 (?)
end

text(.0022,.562,'\textbf{I}')
text(.0067,.56,'\textbf{II}')
text(.0215,.562,'\textbf{III}')

helppoint(1) = min(pib.d(pib.stability(:,1)==1)); %mssing point leftmost on top due to num. inac.
helppoint(2) = Luu(5)-(Luu(6)-Luu(5))*6;

ymax = 800;
ax2 = axes('Position',[.54 .56 .34 .3]);
box on
hold on
patch([pib.d(pib.stability(:,1)==1),fliplr(pib.d(upper)),helppoint(1)],[Lul;flipud(Luu(5:end));helppoint(2)],C.bireg,'facealpha',.2,'edgealpha',0);
plot([0 max(pib.d)],[0 0],'-','linewidth',1.5,'color',C.stable);
plot(bis.d(pib.stability(:,1)==1),Lul,'-','linewidth',1,'color',C.unstable)
plot([helppoint(1) bis.d(upper)],[helppoint(2);Luu(5:end)],'-','linewidth',1,'color',C.unstable)
%plot([pfb.D(1),min(bis.d(upper))],[Lul(1) Luu(5)],'-','linewidth',1,'color',C.unstable)
plot(bis.d(pib.stability(:,2)==-1),Ls,'-','linewidth',1,'color',C.stable,'linewidth',1.5)

text(.05,ymax,'Learning rate $(L_{0\rightarrow 1})$','HorizontalAlignment','center','Verticalalignment','bottom')
% for i = 1:4
%     plot([1 1]*(i)*.02,[0 4000],'color',[.15 .15 .15 .15])
% end
for i = 1:10
    plot(bis.d(pib.stability(:,1)==1),Lgrid(:,i),'--','color',[.15 .15 .15 .5])
end
patch([0 pfb.D(1) pfb.D(1) 0],[0 0 ymax ymax],[.15 .15 .15],'edgealpha',0,'facealpha',.15)
ylim([0 ymax])

text(.092,75,'\textbf{I}','horizontalalignment','center')
text(.092,250,'\textbf{II}','horizontalalignment','center')
text(.092,450,'\textbf{III}','horizontalalignment','center')


set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% Population given learning rate'
figure
names{get(gcf,'number')} = 'Popeq_optstrat';

carry_d = (pfb.data{3}-pfb.data{5})/pfb.data{1};
carry_ay = (pfbay.data{2}-pfbay.data{5})/pfbay.data{1};
carry_mu = (pfbmu.data{3}-[0,1])/pfbmu.data{1};

Bi_d = pfb.D(pfb.xc>0);
sad_d = pfb.U(pfb.xc>0,:);
xc_d = pfb.xc(pfb.xc>0);

ploss_d = sad_d./pfb.X(pfb.xc>0);

Bi_ay = pfbay.Ay(pfbay.xc>0);
sad_ay = pfbay.U(pfbay.xc>0,:);
xc_ay = pfbay.xc(pfbay.xc>0);

ploss_ay = sad_ay./pfbay.X(pfbay.xc>0);

Bi_mu = pfbmu.Mu(pfbmu.xc>0);
sad_mu = pfbmu.U(pfbmu.xc>0,:);
xc_mu = pfbmu.xc(pfbmu.xc>0);

ploss_mu = sad_mu./pfbmu.X(pfbmu.xc>0);

ymax = 1000;
ymaxsoc = .4;

tiledlayout(3,3, 'TileIndexing', 'columnmajor')
nexttile([2 1])
patch([Bi_d(1) Bi_d(end) Bi_d(end) Bi_d(1)],[-ymax -ymax 2*ymax 2*ymax],C.bireg,'facealpha',.2,'edgealpha',0);
hold on
plot([Bi_d(1) Bi_d(1)],[-ymax 2*ymax],'color',C.bireg)
plot([Bi_d(end) Bi_d(end)],[-ymax ymax*.85],'color',C.bireg)

plot(pfb.D,pfb.X(:,1),'color',C.nonmig,'linewidth',1.5)
plot(pfb.D,pfb.X(:,2),'color',C.mig,'linewidth',1.5)
plot(Bi_d,sad_d(:,1),'--','color',C.nonmig,'linewidth',1.5)
plot(Bi_d,sad_d(:,2),'--','color',C.mig,'linewidth',1.5)
plot(Bi_d,xc_d,'-','color',C.nonmig,'linewidth',1.5)
plot([0 pfb.D(end)],[carry_d,carry_d],'k','linewidth',1.5)
%plot(.02:.004:.028,repmat(1000,1,3),'^','color',C.mig,'markerfacecolor',C.mig)
%plot([pfbay.data{4} pfbay.data{4}],[0,ymax],'k')
text(.03,ymax,[strcat('$d_1$ = ',num2str(pfb.data{2}));strcat("$\mu$ = ",num2str(pfb.data{5}))],'HorizontalAlignment', 'right','Verticalalignment','top');
ylim([0 ymax])
xlim([0 .03])
xlabel('Learning rate coefficient $(c)$')
ylabel('Population equilibrium')


nexttile
plot([0 Bi_d(1) Bi_d Bi_d(end) .03],[0;0;1-ploss_d(:,1);1;1],'color',C.nonmig,'linewidth',1.5)
hold on
plot([0 Bi_d(1) Bi_d Bi_d(end) .03],[0;0;1-ploss_d(:,2);1;1],'color',C.mig,'linewidth',1.5)
ylim([0 1])
xlim([0 .03])
xlabel('Learning rate coefficient $(c)$')
ylabel(["Critical population";" fraction"])
set(gca,'ColorOrderIndex',1);
colororder(C.stable);
yyaxis right
plot(pib.d,(pib.eq(:,2)-pib.eq(:,1))./pib.eq(:,2),':','linewidth',1.5)
ylim([0 ymaxsoc]);

nexttile([2 1])
ptch = patch([Bi_ay(1) Bi_ay(end) Bi_ay(end) Bi_ay(1)],[-ymax -ymax 2*ymax 2*ymax],C.bireg,'facealpha',.2,'edgealpha',0);
hold on
plot([Bi_ay(1) Bi_ay(1)],[-ymax ymax],'color',C.bireg)
plot([Bi_ay(end) Bi_ay(end)],[-ymax ymax*.85],'color',C.bireg)

pltnm = plot(pfbay.Ay,pfbay.X(:,1),'color',C.nonmig,'linewidth',1.5);
pltm = plot(pfbay.Ay,pfbay.X(:,2),'color',C.mig,'linewidth',1.5);
plot(Bi_ay,sad_ay(:,1),'--','color',C.nonmig,'linewidth',1.5)
plot(Bi_ay,sad_ay(:,2),'--','color',C.mig,'linewidth',1.5)
plot(Bi_ay,xc_ay,'-','color',C.nonmig,'linewidth',1.5)
pltc = plot([0 .001],[carry_ay,carry_ay],'color',[0 0 0],'linewidth',1.5);
%plot(.00025,1000,'^','color',C.mig,'markerfacecolor',C.mig)
%plot([pfb.data{2} pfb.data{2}],[0,ymax],'k')
text(.00085,ymax,[strcat('$c$ = ',num2str(pfbay.data{3}));strcat("$\mu$ = ",num2str(pfbay.data{5}))],'HorizontalAlignment', 'right','Verticalalignment','top');
ylim([0 ymax])
xlim([.0001 .001])
xlabel('Dist. CNDD $(d_1)$')

dummy1 = plot(nan,nan,'k-');
dummy2 = plot(nan,nan,'k--');
nothing = plot(nan,nan,'color','none');
legend([pltnm pltm pltc nothing nothing nothing dummy1 dummy2 ptch],{'Non-migrants ($\sigma^*>0$)','Migrants\hspace{20.5pt} ($\sigma^*>0$)','Non-migrants ($\sigma^*=0$)',...
    '','','','Stable eq.','Saddle point','Bi-stable region'},'numcolumns',3,'location','northoutside')
set(gca,'xtick',[1 3 5 7]*1e-4)
set(gca,'xticklabel',{1 3 5 7})
xlim([.0001 .00085])


nexttile
plot([Bi_ay Bi_ay(end) .001],[1-ploss_ay(:,1);0;0],'color',C.nonmig,'linewidth',1.5)
hold on
plot([Bi_ay Bi_ay(end) .001],[1-ploss_ay(:,2);0;0],'color',C.mig,'linewidth',1.5)
ylim([0 1])
xlim([.0001 .001])
xlabel('Dist. CNDD $(d_1)$')
% lgd = legend({'Non-migrants','Migrants'},'location','northoutside');
% lgd.Position(1) = lgd.Position(1) + .02;
% lgd.Position(2) = lgd.Position(2)*4/5;
set(gca,'ColorOrderIndex',1);
colororder(C.stable);
yyaxis right
plot(pibay.ay,(pibay.eq(:,2)-pibay.eq(:,1))./pibay.eq(:,2),':','linewidth',1.5)
ylim([0 ymaxsoc]);
annotation('textbox',[0.575,0.38,0,0],'string','$\times 10^{-4}$')
set(gca,'xtick',[1 3 5 7]*1e-4)
set(gca,'xticklabel',{1 3 5 7})
xlim([.0001 .00085])

nexttile([2 1])
patch([Bi_mu(1) Bi_mu(end) Bi_mu(end) Bi_mu(1)],[-ymax -ymax 2*ymax 2*ymax],C.bireg,'facealpha',.2,'edgealpha',0);
hold on
plot([Bi_mu(1) Bi_mu(1)],[-ymax 2*ymax],'color',C.bireg)
plot([Bi_mu(end) Bi_mu(end)],[-ymax ymax*.85],'color',C.bireg)

plot(pfbmu.Mu,pfbmu.X(:,1),'color',C.nonmig,'linewidth',1.5)
plot(pfbmu.Mu,pfbmu.X(:,2),'color',C.mig,'linewidth',1.5)
plot(Bi_mu,sad_mu(:,1),'--','color',C.nonmig,'linewidth',1.5)
plot(Bi_mu,sad_mu(:,2),'--','color',C.mig,'linewidth',1.5)
plot(Bi_mu,xc_mu,'-','color',C.nonmig,'linewidth',1.5)
plot([0 1],carry_mu,'color',[0 0 0],'linewidth',1.5)
%plot(.415:.02:.485,repmat(1000,1,4),'^','color',C.mig,'markerfacecolor',C.mig)
%plot([pfb.data{5} pfb.data{5}],[0,ymax],'k')
text(.625,ymax,[strcat('$c$ = ',num2str(pfbmu.data{4}));strcat("$d_1$ = ",num2str(pfbmu.data{2}))],'HorizontalAlignment', 'right','Verticalalignment','top');
ylim([0 ymax])
xlim([.4 .625])
xlabel('Mortality rate $(\mu)$')
annotation('textbox',[0.575,0.098,0,0],'string','$\times 10^{-4}$')



nexttile
plot([.4 Bi_mu(1) Bi_mu Bi_mu(end) 1],[1;1;1-ploss_mu(:,1);0;0],'color',C.nonmig,'linewidth',1.5)
hold on
plot([.4 Bi_mu(1) Bi_mu Bi_mu(end) 1],[1;1;1-ploss_mu(:,2);0;0],'color',C.mig,'linewidth',1.5)
ylim([0 1])
xlim([.4 .625])
xlabel('Mortality rate $(\mu)$')
set(gca,'ColorOrderIndex',1);
colororder(C.stable);
yyaxis right
plot(pibmu.mu,(pibmu.eq(:,2)-pibmu.eq(:,1))./pibmu.eq(:,2),':','linewidth',1.5)
ylim([0 ymaxsoc]);
ylabel('Critical sociality fraction')

set(gcf,'Position',[350 150 1200 800])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% Phase plot of critical values in strat
figure 
names{get(gcf,'number')} = 'CriticalVals';

imagesc(log(pha.D),pha.Mu,log10(pha.Ay))
xlabel('Learning rate coefficient log (log$_{10}(c)$)')
ylabel('Mortality rate ($\mu$)')
set(gca,'ydir','normal')
cb = colorbar;
yl = ylabel(cb,'Dist. CNDD log (log$_{10}(d_1)$)','FontSize',14,'Rotation',270,'interpreter','latex');
yl.Position(1) = yl.Position(1) + 1.5;

hold on
Traces = 10.^-[99 6 5 4 3 2 1];
Perims = cell(length(Traces),1);
for i = 1:length(Traces)
    Perims{i} = bwperimtrace(pha.Ay > Traces(i),log([pha.D(1) pha.D(end)]),[pha.Mu(1) pha.Mu(end)]);
    plot(Perims{i}{1}(:,1),Perims{i}{1}(:,2),'k','linewidth',1.5)
end  

text(-6,.67,'$<10^{-6}$')
text(-5,.6,'$<10^{-5}$')
text(-3.2,.55,'$<10^{-4}$')
text(-1.7,.5,'$<10^{-3}$')
text(0.2,.44,'$<10^{-2}$')
text(2.8,.38,'$<10^{-1}$')
text(5.2,.33,'$<1$')


set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
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

bitau = pp1d.tau((pp1d.TotPopb-pp1d.TotPops)>1);
bipops = pp1d.TotPops((pp1d.TotPopb-pp1d.TotPops)>1);

colororder([C.totpop;C.mig])
hold on
yyaxis left
ymax = pp1d.TotPopb(pp1d.monind)*1.05;
pat = patch([min(bitau) min(bitau) max(bitau) max(bitau)],[ymin-1 ymax*1.1  ymax*1.1 ymin-1],C.bireg);
pat.EdgeColor = C.bireg;
pat.FaceAlpha = 0.2;

plt1 = plot(pp1d.tau,pp1d.TotPopb,'linewidth',1.5);
plt2 = plot(bitau,bipops,'--','linewidth',1.5);
chos = plot(pp1d.pol,pp1d.TotPopb(pp1d.polind),'o','markersize',5,'color',C.stable,'linewidth',2);
opt = plot(pp1d.mon,pp1d.TotPopb(pp1d.monind),'o','markersize',5,'color',C.opt,'linewidth',2);
ylabel('Total population')
xlabel('Sociality ($\sigma$)')
xlim([0 max(pp1d.tau)])
ylim([0 ymax])
yyaxis right
plot(pp1d.tau,pp1d.yPopb./(pp1d.xPopb+pp1d.yPopb),'linewidth',1.5)
plot([bitau(1),bitau(end)],[0 0],'--','linewidth',1.5)
%plot(pp1d.tau,pp1d.yPops./(pp1d.xPops+pp1d.yPops),'--','linewidth',1.5)
ylim([-.005 1])
ylabel('Migrating fraction of population')


%ylabel('Migrating fraction of population')
legend([ plt1 plt2 chos opt pat],{'Population eq.','Lower equilibrium','Chosen sociality','Ideal sociality','Bi-stable region'},'location','se')

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',17)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%% Population for varied ay and tau
figure
names{get(gcf,'number')} = 'd1Varied';

nump1 = length(popt.tau);
nump2 = length(popt.ay);
tol = 1e-5;

surf(popt.ay,popt.tau,popt.TotPop,popt.prop,'edgecolor','none');
collapse = popt.TotPop; collapse(popt.prop>tol) = nan;
hold on
surf(popt.ay,popt.tau,collapse,popt.prop-.01,'edgecolor','none');
view(45,55)
colormap(customcolormap([0 .15 .3 .5 .99 1],{'#a0ffa0',rgb2hex(C.mig),'#208060',rgb2hex(C.nonmig+[-.1 .1 0]),rgb2hex(C.nonmig),rgb2hex(C.nonmig/2)}))
c = colorbar;
ylabel(c,'Migrating fraction of population','interpreter','latex')

spacing1 = 125;  % better grid size
for i = [1:spacing1:nump1,nump1]
    plot3(popt.ay, ones(nump2,1)*popt.tau(i), popt.TotPop(i,:),'-k');
end
spacing2 = 134;  % better grid size
for i = [1:spacing2:nump2,nump2]
    plot3(ones(nump1,1)*popt.ay(i), popt.tau, popt.TotPop(:,i),'-k');
end

[~,idx] = min(abs(popt.ay-3e-4));
plot3(ones(nump1,1)*3e-4, popt.tau, popt.TotPop(:,idx),'-k','linewidth',2);

Totpol = popt.TotPop(sub2ind([nump1 nump2],popt.polind,1:nump2));

pol0pol = popt.pol(popt.pol == 0);
pol0ax = popt.ay(popt.pol == 0);
pol0tot = Totpol(popt.pol == 0);

pol1pol = popt.pol(popt.pol ~= 0);
pol1polind = popt.polind(popt.pol ~= 0);
pol1ax = popt.ay(popt.pol  ~= 0);
pol1tot = Totpol(popt.pol ~= 0);

dashmidpoint = (pol1ax(end)+pol0ax(1))/2;
dashmidind = length(pol1ax);
dashlength = pol1polind(end);

plot3(ones(dashlength,1)*dashmidpoint, popt.tau(1:dashlength), popt.TotPop(1:dashlength,dashmidind),':','color',C.stable,'linewidth',2.5);

plot3(pol0ax,pol0pol,pol0tot,'color',C.stable,'linewidth',2.5);
p1 = plot3(pol1ax,pol1pol,pol1tot+25,'color',C.stable,'linewidth',2.5);
p2 = plot3(popt.ay,popt.mon,popt.TotPop(sub2ind([nump1 nump2],popt.monind,1:nump2))+25,'color',C.opt,'linewidth',2.5);

xlim([popt.ay(1) popt.ay(end)])
ylim([popt.tau(1) popt.tau(end)])
xlab = xlabel('Dist. CNDD ($d_{1}$)');
% xlab.Rotation = -atan(320/190)*180/(2*pi);
% xlab.Position(2) = -.03;
ylab = ylabel('Sociality ($\sigma$)');

%ylab.Rotation = atan(320/190)*180/(2*pi);
%ylab.Position(1) = -0.03;
%ylab.Position(2) = .3;
zlabel('Total population number')

legend([p1,p2],{'Chosen ssociality','Ideal sociality'},'location','ne')
set(gca,'Xtick',[1 3 5 7 9]*1e-4)
set(gca,'XTickLabel',{1 3 5 7 9})
annotation('textbox',[.06 .32 0 0],'string','$\times 10^{-4}$')
set(gca,'Ztick',[0 2500 5000])

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',17)
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

% path = 'C:/Users/thekn/Pictures/Article1/';
% 
% for i = 1:length(names)
%     if ishandle(i)
%         figure(i)
%         export_fig([path,names{i}],'-png','-transparent','-m5')
%     end
% end
