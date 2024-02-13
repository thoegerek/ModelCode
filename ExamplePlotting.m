%%
load('PIP.mat')
load('Bifurcation_tau.mat')
load('PIP123.mat')
load('PiB_d.mat')
load('Pop_tau.mat')
load('AXAY.mat')

load('PiB_dimless.mat')
load('Pha.mat');
load('Hys.mat');

%%
abmin = 0.0467; %found manually
abmax = 0.0546;
Lower = polyshape([abmin abmin abmax abmax],[pip.tau(1) abmin abmax pip.tau(1)]);
Upper = polyshape([abmin abmin abmax abmax],[pip.tau(end) abmin abmax pip.tau(end)]);

pipd0 = 1.25;
%%
figure(1)
imshow((pip.f>0),'xdata',pip.tau,'ydata',pip.tau,'InitialMagnification',10);
set(gca,'visible','on','ydir','normal')
xlabel('\tau_{r}')
ylabel('\tau_{i}')
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(pip.eqs(pip.stability<0),pip.eqs(pip.stability<0),'bo','markersize',10,'linewidth',2);
ueq = plot(pip.eqs(pip.stability>=0),pip.eqs(pip.stability>=0),'ro','markersize',10,'linewidth',2);
legend([dummy2,dummy1,seq,ueq],{'Positive invasion fitness','Negative invasion fitness','Stable equilibrium','Unstable equilibrium'})

%%
figure(2)
p = plot(bif.tux,bif.ux,'.r',bif.tsx,bif.sx,'.k',bif.tuy,bif.uy,'.m',bif.tsy,bif.sy,'.b');
hold on
legend([p(2),p(1),p(4),p(3)],{'Stable x*','Unstable x*','Stable y*','Unstable y*'})
xlim([bif.tau(1) bif.tau(end)])
xlabel('\tau')
ylabel('Populatoin size')
%%
figure(3)
axis;
hold on
for i = 1:size(pib.eq,2)
    plot(pib.d(pib.stability(:,i)==-1),pib.eq(pib.stability(:,i)==-1,i),'b-','linewidth',1.5)
    plot(pib.d(pib.stability(:,i)==1),pib.eq(pib.stability(:,i)==1,i),'r--','linewidth',1.5)
end
abc = {'a','b','c'};
for i = 1:length(ds)
    plot([ds ds],[0 1],':k')
    text(ds(i)+.007,0.993,abc(i))
end
plot([0.125 0.125],[0 1],':','color',[.7 .7 .7])
xlim([0 1])
ylim([0 1])
xlabel('d')
ylabel('\tau')
legend('Stable strategy equlibrium','Unstable strategy equilibrium','Points of reference')
%%
fig = figure(4);
pip123 = {pip1,pip2,pip3};
for i = 1:3
    subplot(1,3,i)
    imshow((pip123{i}.f>0),'xdata',pip123{i}.tau,'ydata',pip123{i}.tau,'InitialMagnification',10);
    set(gca,'visible','on','ydir','normal')
    text(0.1,1.5,abc(i),'color','w')
    xlabel('\tau_{r}')
    if i == 1
        ylabel('\tau_{i}')
    end
    hold on
    dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
    dummy2 = plot(inf,inf,'squarek');
    seq = plot(pip123{i}.eqs(pip123{i}.stability<0),pip123{i}.eqs(pip123{i}.stability<0),'bo','markersize',5,'linewidth',1.5);
    ueq = plot(pip123{i}.eqs(pip123{i}.stability>=0),pip123{i}.eqs(pip123{i}.stability>=0),'ro','markersize',5,'linewidth',1.5);
    if i ==2
        legend([dummy2,dummy1,seq,ueq],{'Positive invasion fitness','Negative invasion fitness','Stable equilibrium','Unstable equilibrium'})
    end
end
%%
nump1 = length(popt.tau);
nump2 = length(popt.ax);
figure(5)
surf(popt.ax,popt.tau,popt.TotPop,popt.prop,'edgecolor','none');
view(45,45)
colormap('jet')
c = colorbar;
ylabel(c,'{y*}/{x*}')
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
xlabel('a_{x}')
ylabel('\tau')
ylim([popt.tau(1) popt.tau(end)])
zlabel('x* + y*')

legend([p1,p2],{'Chosen strategy','Ideal strategy'},'location','ne')
%%
nump1 = length(axay.ax);
nump2 = length(axay.ay);
figure(6)
surf(axay.ax,axay.ay,log(axay.TotPop),log(axay.prop+1),'edgecolor','none');
view(135,45)
colormap('jet')
c = colorbar;
ylabel(c,'ln({y*}/{x*} +1)')
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
xlabel('a_{x}')
ylabel('a_{y}')
ylim([axay.ay(1) axay.ay(round(end/3))])
zlabel('ln(x* + y* +1)')
%%



%%
figure(7)
axis;
hold on
for i = 1:size(pidi.stability,2)
    plot(pidi.d(pidi.stability(:,i)==1),pidi.eq(pidi.stability(:,i)==1,i),'m.')
    plot(pidi.d(pidi.stability(:,i)==-1),pidi.eq(pidi.stability(:,i)==-1,i),'k.')
end
xlim([0 max(pidi.d)])
xlabel('d \cdot b^{3} \cdot c^{-2} \cdot a_{x}^{-1}');
ylabel('\tau [b/c]')
%%
figure(8)
imagesc(log(pha.bif),'xdata',log(pha.Ay),'ydata',pha.Mu)
xlabel('ln(a_{y})')
ylabel('\mu')
set(gca,'ydir','normal')
cb = colorbar;
yl = ylabel(cb,'d \cdot b^{3} \cdot c^{-2} \cdot a_{x}^{-1}','FontSize',12,'Rotation',270);
yl.Position(1) = yl.Position(1) + 1.5;
%%
figure(9)
yyaxis left
plot(1:length(hys.D),hys.D,'k',1:length(hys.D),hys.tau,'--b')
yyaxis right
plot(1:length(hys.D),hys.eq)
legend('D','\tau*','x','y')