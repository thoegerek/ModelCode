%%
close all
%%
load('PIP.mat')

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
%%
figure
names{get(gcf,'number')} =  'PIP';

tiledlayout(3,3)
nexttile([3 3])
imshow((pip.f(:,:,2)<0),'xdata',pip.tau,'ydata',pip.tau)
set(gca,'visible','on','ydir','normal','box','off')
ylabel('Invader sociality ($\sigma_{m}$)','fontsize',12)
xlabel('Resident sociality ($\sigma_{n}$)','fontsize',12)
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(0,0,'o',pip.eqs(end),pip.eqs(end),'o','color',C.stable,'markersize',10,'linewidth',2);
%ueq = plot(pip.eqs(pip.stability>=0),pip.eqs(pip.stability>=0),'-o','color',C.unstable,'markersize',10,'linewidth',2);


set(0,'defaulttextInterpreter','tex') 

text(.40,.1,'$-$','color','k')

text(.33,.5,'$+$','color','w')
text(.1095,.07,'$-$','fontsize',16)

text(.185,.6,'$-$','fontsize',24)
text(.185,.1,'$+$','color','w','fontsize',24)
text(.8,.3,'$+$','color','w','fontsize',24)



plot(1,.8,'^','color',C.bireg,'markerfacecolor',C.bireg,'markersize',8);
plot(.8,.6,'s','color',C.bireg,'markerfacecolor',C.bireg,'markersize',10);
plot(.6,.4,'diamond','color',C.bireg,'markerfacecolor',C.bireg,'markersize',9);
plot(.35,.2,'^','color',[0 .5 .1],'markerfacecolor',[0 .5 .1],'markersize',8);
plot(.6,.7,'s','color',[0 .5 .1],'markerfacecolor',[0 .5 .1],'markersize',10);

legend([seq(1),dummy1,dummy2],{'Evolutionary Stable Strategy','Positive invasion fitness','Negative invasion fitness'...
    },'location','ne','fontsize',12);

set(gcf,'Position',[680   418   713   680])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

export_fig('C:/Users/thekn/Pictures/Patritation/PIP','-png','-transparent','-m5')