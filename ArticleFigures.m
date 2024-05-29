%%

load('Invasion.mat')
load('PIP.mat')
load('PIP2.mat')
load('PIP3.mat')
load('PIP4.mat')
load('Bifurcation_tau.mat')
load('PhasePop.mat')
load('PiB_d.mat')
load('PopFromBif.mat')
load('PopFromBif_ay.mat')
load('PopFromBif_mu.mat')
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
C.opt = [0.85 0.2 0.85];
C.totpop = [0 0 0];
C.nonmig = [0.65 0.45 0];
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
% styles = flip({'^','^','^','^','^','^'...
%     's','s','s','s',...
%     'diamond','diamond','diamond','diamond',...
%     'o','o'});

styles = flip({'-','-.','-.','-.','-.','-.','-.',...
'-','-','-','-',...
'-.','-.','-.','-.',...
'-','-','-','-','-','-','-','-','-','-'});

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

plot(inv.T2(find(max(inv.X2(:,3:4),[],2)>=ymin,1,'first')),ymin,'^','color',C.bireg,'markerfacecolor',C.bireg,'markersize',6);
plot(inv.T3(find(max(inv.X3(:,3:4),[],2)>=ymin,1,'first')),ymin,'s','color',C.bireg,'markerfacecolor',C.bireg,'markersize',8);
plot(inv.T4(find(max(inv.X4(:,3:4),[],2)>=ymin,1,'first')),ymin,'diamond','color',C.bireg,'markerfacecolor',C.bireg,'markersize',7);


text(25,max(inv.X1,[],'all')+40,'$\sigma = 1$')
text(550,max(inv.X2,[],'all')+40,'$\sigma = 0.8$')
text(1250,max(inv.X4,[],'all')+40,'$\sigma = 0.6$')
text(2200,max(inv.X5,[],'all')+40,'$\sigma = 0.4$')
text(2232,1850,'$\mu = 5.7 -$')

patch([inv.T5(1) inv.T5(1) inv.T5(end) inv.T5(end)],[0 ymax ymax 0],[.5 .5 .5],'facealpha',.3,'edgealpha',0)

xlim([0,Tcat(end)])
ylim([ymin ymax])
legend([plt(8),plt(9),plt(1)],{'Non-migrating population','Migrating population','Total population',...
    },'location','nw');
xlabel('Time (years)')
ylabel('Population number')


set(gcf,'Position',[550 350 825 575])
set(findall(gcf,'-property','FontSize'),'FontSize',14)
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

text(.2,.6,'$-$','fontsize',24)
text(.2,.1,'$+$','color','w','fontsize',24)
text(.8,.3,'$+$','color','w','fontsize',24)

plot([.4 .6 .6 .8 .8 1],[.4 .4 .6 .6 .8 .8],'--','color',[.5 .5 .5],'linewidth',1.5);
plot(.5,.4,'<','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])
plot(.6,.5,'v','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])
plot(.7,.6,'<','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])
plot(.8,.7,'v','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])
plot(.9,.8,'<','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5])

example = plot(nan,nan,'<--','color',[.5 .5 .5],'linewidth',1.5,'markerfacecolor',[.5 .5 .5]);

plot(1,.8,'^','color',C.bireg,'markerfacecolor',C.bireg,'markersize',6);
plot(.8,.6,'s','color',C.bireg,'markerfacecolor',C.bireg,'markersize',8);
plot(.6,.4,'diamond','color',C.bireg,'markerfacecolor',C.bireg,'markersize',7);

legend([seq(1),example],{'Evolutionary Stable Strategy','Path of adaptation (fig. 1)'...
    },'location','ne','fontsize',12);


nexttile
imshow(((pip4.f(:,:,1)<0)+(pip4.f(:,:,2)<0)*2)/3,'xdata',pip4.tau,'ydata',pip4.tau)
hold on
set(gca,'visible','on','ydir','normal','box','off')
plot(0,0,'o','color',C.stable,'markersize',7,'linewidth',1.5)
text(.025,.9,'$d_1\uparrow$','color','w')
%
nexttile
imshow(((pip2.f(:,:,1)<0)+(pip2.f(:,:,2)<0)*2)/3,'xdata',pip2.tau,'ydata',pip2.tau)
hold on
set(gca,'visible','on','ydir','normal','box','off')
plot(0,0,'o',max(pip2.eqs),max(pip2.eqs),'o','color',C.stable,'markersize',7,'linewidth',1.5)
xlabel('Resident sociality ($\sigma_{n}$)','fontsize',12)
text(.025,.9,'$d_1\uparrow$','color','w')
text(.025,.75,'$c\uparrow$','color','w')
%
nexttile
imshow(((pip3.f(:,:,1)<0)+(pip3.f(:,:,2)<0)*2)/3,'xdata',pip3.tau,'ydata',pip3.tau)
hold on
set(gca,'visible','on','ydir','normal','box','off')
plot(0,0,'o',max(pip3.eqs),max(pip3.eqs),'o','color',C.stable,'markersize',7,'linewidth',1.5)

text(.025,.9,'$c\uparrow$','color','w')
text(.025,.75,'$\lambda\downarrow$','color','w')


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

text(.02,1100,'(a)')
text(.355,1100,'(b)')
text(.935,1100,'(c)')

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

sols = plot(nan,nan,'color',[.5 .5 .5]);%plot(A1(:,1),A1(:,2),'-k');
manif = plot(nan,nan,'k','linewidth',1.5);
%plot(A2(:,1),A2(:,2),'-k')
%plot(A3(:,1),A3(:,2),'-k')

quiver(php.X,php.Y,php.dx,php.dy,'linewidth',1.5, 'AutoScaleFactor',0.5,'color',[.8 .7 .2])

psta = plot(php.steqx,php.steqy,'o','color',C.stable,'markersize',10,'linewidth',2);
pusta = plot(php.usteqx(1),php.usteqy(1),'o','color',C.unstable,'markersize',10,'linewidth',2);
psad = plot(php.usteqx(2),php.usteqy(2),'o','color',C.opt,'markersize',10,'linewidth',2);
xlabel('Non-migrant population ($n_{0}$)')
ylabel('Migrant population ($n_{1}$)')
legend([psta,psad,pusta,sols,manif,red,green],{'Stable eq.','Saddle point','Unstable eq.','Random trajectories','Attracting manifold'...
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
ymax = .8;
axis;
hold on
patch([pib.d(pib.stability(:,1)==1),fliplr(pib.d(upper))],[pib.eq(pib.stability(:,1)==1,1);flipud(tauupper)],C.bireg,'facealpha',.2,'edgealpha',0);
ptch = patch(nan,nan,C.bireg,'edgecolor',C.unstable,'facealpha',.2);
plt1 = [];
plt2 = [];
plot([0 max(pib.d)],[0 0],'-','linewidth',1.5,'color',C.stable);
plot(bis.d(upper),tauupper,'-','linewidth',1,'color',C.unstable)
plot([bis.d(upper(1)-4),bis.d(upper(1))],[max(pib.eq(pib.stability(:,i)==1,i)) tauupper(1)],'-','linewidth',1,'color',C.unstable)
for i = size(pib.eq,2):-1:1
    plt1{i} = plot(pib.d(pib.stability(:,i)==-1),pib.eq(pib.stability(:,i)==-1,i),'-','linewidth',1.5,'color',C.stable);
    plt2{i} = plot(pib.d(pib.stability(:,i)==1),pib.eq(pib.stability(:,i)==1,i),'-','linewidth',1,'color',C.unstable);
end


plot(1,0,'o','color',C.stable,'markersize',5,'linewidth',2);

%plot([pips.d pips.d],[ymin ymax],':','color',[.7 .7 .7])
%text(pips.d+.005,.04,'-figure 5 (top-left)','color',[.5 .5 .5])

plot([1 1],[ymin ymax],':','color',[.5 .5 .5])
text(1+.005,.03,'-figure 3','color',[.3 .3 .3])

xlim([0 .1])
ylim([ymin ymax])
xlabel('Learning rate coefficient (c)')

ylabel('Sociality ($\sigma$)')


% xzoom = [.005 .025];
% yzoom = [.3 .5];
% rectangle('Position', [xzoom(1), yzoom(1), xzoom(2)-xzoom(1),  yzoom(2)-yzoom(1)]);
% plot([xzoom(2) .08],[.45 .45],'k-')
legend([plt1{2},ptch],{'Stable strategy equlibrium','Bistable region'},'location','ne')
% ax1 = gca;

% ax2 = axes('Position',[.6 .6 .3 .3]);
%     ptch = patch([pib.d(pib.stability(:,1)==1),fliplr(pib.d(upper))],[pib.eq(pib.stability(:,1)==1,1);flipud(tauupper)],C.bireg,'facealpha',.2,'edgealpha',0);
%     hold on
%     for i = size(pib.eq,2):-1:1
%         plot(pib.d(pib.stability(:,i)==-1),pib.eq(pib.stability(:,i)==-1,i),'-','linewidth',1.5,'color',C.stable);
%         plot(pib.d(pib.stability(:,i)==1),pib.eq(pib.stability(:,i)==1,i),'--','linewidth',1.5,'color',C.unstable);
%     end
%     plot([0 max(pib.d)],[0 0],'-','linewidth',1.5,'color',C.stable);
%     plot(bis.d(upper),tauupper,'--','linewidth',1.5,'color',C.unstable)
%     plot([1 1],[ymin ymax],':','color',[.7 .7 .7])
%     plot(1,pip.eqs(pip.stability<0),'o','color',C.stable,'markersize',5,'linewidth',2);
%     plot([1 1],pip.eqs(pip.stability>=0),'o-','color',C.unstable,'markersize',5,'linewidth',2);
%     rectangle('Position', [xzoom(1), yzoom(1), xzoom(2)-xzoom(1),  yzoom(2)-yzoom(1)]);
%     xlim(xzoom)
%     ylim(yzoom)






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

tiledlayout(3,3, 'TileIndexing', 'columnmajor')
nexttile([2 1])
patch([Bi_d(1) Bi_d(end) Bi_d(end) Bi_d(1)],[-ymax -ymax 2*ymax 2*ymax],C.bireg,'facealpha',.2,'edgecolor',C.bireg);
hold on
plot(pfb.D,pfb.X(:,1),'color',C.nonmig,'linewidth',1.5)
plot(pfb.D,pfb.X(:,2),'color',C.mig,'linewidth',1.5)
plot(Bi_d,sad_d(:,1),'--','color',C.nonmig,'linewidth',1.5)
plot(Bi_d,sad_d(:,2),'--','color',C.mig,'linewidth',1.5)
plot(Bi_d,xc_d,'-','color',C.nonmig,'linewidth',1.5)
plot([0 pfb.D(end)],[carry_d,carry_d],'color',C.stable,'linewidth',1.5)
plot(.02:.004:.028,repmat(1000,1,3),'^','color',C.mig,'markerfacecolor',C.mig)
%plot([pfbay.data{4} pfbay.data{4}],[0,ymax],'k')
text(.03,975,[strcat('$d_1$ = ',num2str(pfb.data{2}));strcat("$\mu$ = ",num2str(pfb.data{5}))],'HorizontalAlignment', 'right','Verticalalignment','top');
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

nexttile([2 1])
ptch = patch([Bi_ay(1) Bi_ay(end) Bi_ay(end) Bi_ay(1)],[-ymax -ymax 2*ymax 2*ymax],C.bireg,'facealpha',.2,'edgecolor',C.bireg);
hold on
pltnm = plot(pfbay.Ay,pfbay.X(:,1),'color',C.nonmig,'linewidth',1.5);
pltm = plot(pfbay.Ay,pfbay.X(:,2),'color',C.mig,'linewidth',1.5);
plot(Bi_ay,sad_ay(:,1),'--','color',C.nonmig,'linewidth',1.5)
plot(Bi_ay,sad_ay(:,2),'--','color',C.mig,'linewidth',1.5)
plot(Bi_ay,xc_ay,'-','color',C.nonmig,'linewidth',1.5)
pltc = plot([0 .001],[carry_ay,carry_ay],'color',C.stable,'linewidth',1.5);
plot(.00025,1000,'^','color',C.mig,'markerfacecolor',C.mig)
%plot([pfb.data{2} pfb.data{2}],[0,ymax],'k')
text(.001,975,[strcat('$c$ = ',num2str(pfbay.data{3}));strcat("$\mu$ = ",num2str(pfbay.data{5}))],'HorizontalAlignment', 'right','Verticalalignment','top');
ylim([0 ymax])
xlim([.0002 .001])
xlabel('Dist. CNDD $(d_1)$')

dummy1 = plot(nan,nan,'k-');
dummy2 = plot(nan,nan,'k--');
nothing = plot(nan,nan,'color','none');
legend([pltnm pltm pltc nothing nothing nothing dummy1 dummy2 ptch],{'Non-migrants ($\sigma*>0$)','Migrants\hspace{20.5pt} ($\sigma*>0$)','Non-migrants ($\sigma*=0$)',...
    '','','','Stable eq.','Saddle point','Bi-stable region'},'numcolumns',3,'location','northoutside')


nexttile
plot([Bi_ay Bi_ay(end) .001],[1-ploss_ay(:,1);0;0],'color',C.nonmig,'linewidth',1.5)
hold on
plot([Bi_ay Bi_ay(end) .001],[1-ploss_ay(:,2);0;0],'color',C.mig,'linewidth',1.5)
ylim([0 1])
xlim([.0002 .001])
xlabel('Dist. CNDD $(d_1)$')
lgd = legend({'Non-migrants','Migrants'},'location','northoutside');
lgd.Position(1) = lgd.Position(1) + .02;
lgd.Position(2) = lgd.Position(2)*4/5;

nexttile([2 1])
patch([Bi_mu(1) Bi_mu(end) Bi_mu(end) Bi_mu(1)],[-ymax -ymax 2*ymax 2*ymax],C.bireg,'facealpha',.2,'edgecolor',C.bireg);
hold on
plot(pfbmu.Mu,pfbmu.X(:,1),'color',C.nonmig,'linewidth',1.5)
plot(pfbmu.Mu,pfbmu.X(:,2),'color',C.mig,'linewidth',1.5)
plot(Bi_mu,sad_mu(:,1),'--','color',C.nonmig,'linewidth',1.5)
plot(Bi_mu,sad_mu(:,2),'--','color',C.mig,'linewidth',1.5)
plot(Bi_mu,xc_mu,'-','color',C.nonmig,'linewidth',1.5)
plot([0 1],carry_mu,'color',C.stable,'linewidth',1.5)
plot(.415:.02:.485,repmat(1000,1,4),'^','color',C.mig,'markerfacecolor',C.mig)
%plot([pfb.data{5} pfb.data{5}],[0,ymax],'k')
text(.625,975,[strcat('$c$ = ',num2str(pfbmu.data{4}));strcat("$d_1$ = ",num2str(pfbmu.data{2}))],'HorizontalAlignment', 'right','Verticalalignment','top');
ylim([0 ymax])
xlim([.4 .625])
xlabel('Mortality rate $(\mu)$')


nexttile
plot([.4 Bi_mu(1) Bi_mu Bi_mu(end) 1],[1;1;1-ploss_mu(:,1);0;0],'color',C.nonmig,'linewidth',1.5)
hold on
plot([.4 Bi_mu(1) Bi_mu Bi_mu(end) 1],[1;1;1-ploss_mu(:,2);0;0],'color',C.mig,'linewidth',1.5)
ylim([0 1])
xlim([.4 .625])
xlabel('Mortality rate $(\mu)$')

set(gcf,'Position',[350 150 1200 800])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
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
legend([ plt1 plt2 chos opt pat],{'Population eq.','Non-migratory eq.','Chosen sociality','Ideal sociality','Bi-stable region'},'location','se')

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
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
colormap(customcolormap([0 .2 .4 .99 1],{'#a0ffa0',rgb2hex(C.mig),rgb2hex((C.nonmig*2+C.mig)/3),rgb2hex(C.nonmig),rgb2hex(C.nonmig/2)}))
c = colorbar;
ylabel(c,'Migrating fraction of population','interpreter','latex')

spacing1 = 125;  % better grid size
for i = [1:spacing1:nump1,nump1]
    plot3(popt.ay, ones(nump2,1)*popt.tau(i), popt.TotPop(i,:),'-k');
end
spacing2 = 125;  % better grid size
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

legend([p1,p2],{'Chosen strategy\,','Ideal strategy'},'location','ne')

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
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
