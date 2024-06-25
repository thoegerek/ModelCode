%%
close all
%%

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
ax = .001;
ay = .0003;
b = 1;
c = .2;
d = .015;
mu = .5;
mubad = .570;

tau1 = 1;
tau2 = .8;
tau3 = .6;
tau4 = .4;
%%
X10 = [0;0;1;1];
[T1,X1] = ode15s(@myModelInvader,[0 100],X10,[],ax,ay,b,c,d,mu,0,tau1);
%%
figure

ymax = 2000;
ymin = 5;
plot(T1,X1,'linewidth',1.5);
colororder([
    repmat([C.nonmig;C.mig],12,1)]);

xlim([0,T1(end)])
ylim([ymin ymax])
 legend('Non-migrating population','Migrating population','location','ne');
xlabel('Time')
ylabel('Population number')


set(gcf,'Position',[550 350 825 575])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

%export_fig('C:/Users/thekn/Pictures/Patritation/sigma1','-png','-transparent','-m5')
X1eq = X1(end,3:4);
T1eq = T1(end);
%%
X20 = [X1eq';1;0];
[T2,X2] = ode15s(@myModelInvader,[0 500],X20,[],ax,ay,b,c,d,mu,tau1,tau2);
%%
figure;
names{get(gcf,'number')} = 'PopInv';


ymax = 2000;
ymin = 5;
plot(T2,X2(:,1:2),'linewidth',1.5);
hold on
plot(T2,X2(:,3:4),'-.','linewidth',1.5);
colororder([
    repmat([C.nonmig;C.mig],12,1)]);

xlim([0,T2(end)])
ylim([ymin ymax])
legend('Non-migrating residents','Migrating residents','Non-migrating invaders','Migrating invaders','location','ne','numcolumns',2);
xlabel('Time')
ylabel('Population number')


set(gcf,'Position',[550 350 825 575])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

%export_fig('C:/Users/thekn/Pictures/Patritation/sigma08','-png','-transparent','-m5')
X2eq = X2(end,3:4);
T2eq = T2(end);
%%
X30 = [X2eq';1;0];
[T3,X3] = ode15s(@myModelInvader,[0 500],X30,[],ax,ay,b,c,d,mu,tau2,tau3);
%%
figure;
names{get(gcf,'number')} = 'PopInv';


ymax = 2000;
ymin = 5;
plot(T3,X3(:,1:2),'-.','linewidth',1.5);
hold on
plot(T3,X3(:,3:4),'linewidth',1.5);
colororder([
    repmat([C.nonmig;C.mig],12,1)]);

xlim([0,T3(end)])
ylim([ymin ymax])
legend('Non-migrating residents','Migrating residents','Non-migrating invaders','Migrating invaders','location','ne','numcolumns',2);
xlabel('Time')
ylabel('Population number')


set(gcf,'Position',[550 350 825 575])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

%export_fig('C:/Users/thekn/Pictures/Patritation/sigma06','-png','-transparent','-m5')
X3eq = X3(end,3:4);
T3eq = T3(end);
%%
X40 = [X3eq';1;0];
[T4,X4] = ode15s(@myModelInvader,[0 1000],X40,[],ax,ay,b,c,d,mu,tau3,tau4);
%%
fig = figure;
names{get(gcf,'number')} = 'PopInv';


ymax = 2000;
ymin = 5;
plot(T4,X4(:,1:2),'linewidth',1.5);
hold on
plot(T4,X4(:,3:4),'-.','linewidth',1.5);
colororder([
    repmat([C.nonmig;C.mig],12,1)]);

xlim([0,T4(end)])
ylim([ymin ymax])
legend('Non-migrating residents','Migrating residents','Non-migrating invaders','Migrating invaders','location','ne','numcolumns',2);
xlabel('Time')
ylabel('Population number')


set(gcf,'Position',[550 350 825 575])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

export_fig('C:/Users/thekn/Pictures/Patritation/sigma04','-png','-transparent','-m5')
X4eq = X4(end,3:4);
T4eq = T4(end);
%%
tau5 = 0.35;
X50 = [0;0;1000;1000];
[~,X5] = ode15s(@myModelInvader,[0 1000],X50,[],ax,ay,b,c,d,mu,0,tau5);
%%
tau6 = 0.2;
X60 = [X5(end,3:4)';100;100];
[T6,X6] = ode15s(@myModelInvader,[0 100],X60,[],ax,ay,b,c,d,mu,tau5,tau6);
%%
fig = figure;
names{get(gcf,'number')} = 'PopInv';


ymax = 2000;
ymin = 5;
plot([0;10+T6],[X5(end,3:4); X6(:,1:2)],'linewidth',1.5);
hold on
plot(10+T6,X6(:,3:4),'-.','linewidth',1.5);
plot(10,X6(1,3),'x','linewidth',2,'markersize',20)
colororder([
    repmat([C.nonmig;C.mig],12,1)]);

xlim([0,T6(end)])
ylim([ymin ymax])
legend('Non-migrating residents','Migrating residents','Non-migrating invaders','Migrating invaders','location','ne','numcolumns',2);
xlabel('Time')
ylabel('Population number')


set(gcf,'Position',[550 350 825 575])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

export_fig('C:/Users/thekn/Pictures/Patritation/sigmaNothing02','-png','-transparent','-m5')
%%
tau7 = 0.6;
X70 = [0;0;1000;1000];
[~,X7] = ode15s(@myModelInvader,[0 1000],X70,[],ax,ay,b,c,d,mu,0,tau7);
%%
tau8 = 0.7;
X60 = [X7(end,3:4)';100;100];
[T8,X8] = ode15s(@myModelInvader,[0 200],X60,[],ax,ay,b,c,d,mu,tau7,tau8);
%%
fig = figure;
names{get(gcf,'number')} = 'PopInv';


ymax = 2000;
ymin = 5;
plot([0;10+T8],[X7(end,3:4); X8(:,1:2)],'linewidth',1.5);
hold on
plot(10+T8,X8(:,3:4),'-.','linewidth',1.5);
plot(10,X8(1,3),'x','linewidth',2,'markersize',20)
colororder([
    repmat([C.nonmig;C.mig],12,1)]);

xlim([0,T8(end)])
ylim([ymin ymax])
legend('Non-migrating residents','Migrating residents','Non-migrating invaders','Migrating invaders','location','ne','numcolumns',2);
xlabel('Time')
ylabel('Population number')


set(gcf,'Position',[550 350 825 575])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

export_fig('C:/Users/thekn/Pictures/Patritation/sigmaNothing07','-png','-transparent','-m5')


