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

tau = .35;
T = (1:.025:5).^2;
%%
X10 = [1500;500];
[~,X1] = ode15s(@myModel,[T 1000],X10,[],ax,ay,b,c,d,mu,tau);
%%
close all
figure(1)
for t = 1:length(T)
    tiledlayout(1,2)
    nexttile
    plot(T(1:t),X1(1:t,1),'color',C.nonmig);
    hold on
    plot(T(1:t),X1(1:t,2),'color',C.mig);
    plot(T(1),X1(1,1),'x','linewidth',1.5,'color',C.nonmig)
    plot(T(1),X1(1,2),'x','linewidth',1.5,'color',C.mig)
    xlim([-.5 T(end)])
    ylim([0 2000])
    legend('Non-migrating population','Migrating population','location','ne');
    ylabel('Population number')
    xlabel('Time')
    nexttile
    plot(X1(1:t,1),X1(1:t,2),'linewidth',1.5,'color',[0 .5 .1])
    hold on
    plot(X1(1,1),X1(1,2),'x','linewidth',1.5,'color',[0 .5 .1])
    xlim([0 2000])
    ylim([0 2000])
    ylabel('Migrants')
    xlabel('Non-migrants')

    set(gcf,'Position',[350 150 1000 500])
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(findall(gcf,'-property','interpreter'),'interpreter','latex')
    
    %exportgraphics(gcf,'C:/Users/thekn/Pictures/Patritation/Pha1.gif',"Append",t~=1)
    %pause(.1)
end
plot([X1(end-1,1) X1(end,1)],[X1(end-1,2) X1(end,2)],'linewidth',1.5,'color',[0 .5 .1])
plot(X1(end,1),X1(end,2),'o','linewidth',1.5,'color',[0 .5 .1])
exportgraphics(gcf,'C:/Users/thekn/Pictures/Patritation/Pha1.gif',"Append",true)
%%
X20 = [1900;1200];
[~,X2] = ode15s(@myModel,[T 1000],X20,[],ax,ay,b,c,d,mu,tau);
upper = X2(end,:);
%%
close all
figure(2)
for t = 1:length(T)
    tiledlayout(1,2)
    nexttile
    plot(T(1:t),X2(1:t,1),'color',C.nonmig);
    hold on
    plot(T(1:t),X2(1:t,2),'color',C.mig);
    plot(T(1),X2(1,1),'x','linewidth',1.5,'color',C.nonmig)
    plot(T(1),X2(1,2),'x','linewidth',1.5,'color',C.mig)
    xlim([-.5 T(end)])
    ylim([0 2000])
    legend('Non-migrating population','Migrating population','location','ne');
    ylabel('Population number')
    xlabel('Time')
    nexttile
    hold on
    plot(X1(1,1),X1(1,2),'x','linewidth',1.5,'color',[.3 .8 .4])
    plot(X1(:,1),X1(:,2),'linewidth',1.5,'color',[0 .5 .1 .5])
    plot(X1(end,1),X1(end,2),'o','linewidth',1.5,'color',[0 .5 .1])

    plot(X2(1,1),X2(1,2),'x','linewidth',1.5,'color',[0 .5 .1])
    plot(X2(1:t,1),X2(1:t,2),'linewidth',1.5,'color',[0 .5 .1])
    xlim([0 2000])
    ylim([0 2000])
    ylabel('Migrants')
    xlabel('Non-migrants')

    set(gcf,'Position',[350 150 1000 500])
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(findall(gcf,'-property','interpreter'),'interpreter','latex')

    %exportgraphics(gcf,'C:/Users/thekn/Pictures/Patritation/Pha2.gif',"Append",t~=1)
    %pause(.1)
end
plot([X2(end-1,1) X2(end,1)],[X2(end-1,2) X2(end,2)],'linewidth',1.5,'color',[0 .5 .1])
exportgraphics(gcf,'C:/Users/thekn/Pictures/Patritation/Pha2.gif',"Append",true)
%%
X30 = [750;50];
[~,X3] = ode15s(@myModel,[T 1000],X30,[],ax,ay,b,c,d,mu,tau);
lower = X3(end,:);
%%
close all
figure(3)
for t = 1:length(T)
    tiledlayout(1,2)
    nexttile
    plot(T(1:t),X3(1:t,1),'color',C.nonmig);
    hold on
    plot(T(1:t),X3(1:t,2),'color',C.mig);
    plot(T(1),X3(1,1),'x','linewidth',1.5,'color',C.nonmig)
    plot(T(1),X3(1,2),'x','linewidth',1.5,'color',C.mig)
    xlim([-.5 T(end)])
    ylim([0 2000])
    legend('Non-migrating population','Migrating population','location','ne');
    ylabel('Population number')
    xlabel('Time')
    nexttile
    hold on
    plot(X1(1,1),X1(1,2),'x','linewidth',1.5,'color',[.3 .8 .4])
    plot(X2(1,1),X2(1,2),'x','linewidth',1.5,'color',[.3 .8 .4])
    plot(X1(end,1),X1(end,2),'o','linewidth',1.5,'color',[0 .5 .1 .5])
    plot(X2(:,1),X2(:,2),'linewidth',1.5,'color',[0 .5 .1 .5])
    plot(X1(:,1),X1(:,2),'linewidth',1.5,'color',[0 .5 .1 .5])

    plot(X3(1,1),X3(1,2),'x','linewidth',1.5,'color',[.8 0 0])
    plot(X3(1:t,1),X3(1:t,2),'linewidth',1.5,'color',[.8 0 0])
    xlim([0 2000])
    ylim([0 2000])
    ylabel('Migrants')
    xlabel('Non-migrants')

    set(gcf,'Position',[350 150 1000 500])
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(findall(gcf,'-property','interpreter'),'interpreter','latex')

    %exportgraphics(gcf,'C:/Users/thekn/Pictures/Patritation/Pha3.gif',"Append",t~=1)
    %pause(.1)
end
plot([X3(end-1,1) X3(end,1)],[X3(end-1,2) X3(end,2)],'linewidth',1.5,'color',[.8 0 0])
plot(X3(end,1),X3(end,2),'o','linewidth',1.5,'color',[.8 0 0])
exportgraphics(gcf,'C:/Users/thekn/Pictures/Patritation/Pha3.gif',"Append",true)
%%
T4 = (1:.025:5).^4;
N = 1000;
Xs0 = rand(N,2)*2000;
Xs = cell(N,1);
Carr = ones(N,1);
for i = 1:N
    [~,Xs{i}] = ode15s(@myModel,T,Xs0(i,:),[],ax,ay,b,c,d,mu,tau);
    if norm(Xs{i}(end,:)-upper,2) > norm(Xs{i}(end,:)-lower,2)
        Carr(i) = 2;
    end
end
%%
close all
figure(4)
for t = 1:length(T4)
    hold off
    plot(X1(:,1),X1(:,2),'linewidth',1.5,'color',[0 .5 .1 .5])
    hold on
    plot(X2(:,1),X2(:,2),'linewidth',1.5,'color',[0 .5 .1 .5])
    plot(X3(:,1),X3(:,2),'linewidth',1.5,'color',[.8 0 0 .5])
    plot(X2(end,1),X2(end,2),'o','linewidth',1.5,'color',[0 .5 .1])
    plot(X3(end,1),X3(end,2),'o','linewidth',1.5,'color',[.8 0 0])

    for i = 1:N
        if Carr(i) == 1
            color = [.1 .7 .3 .25];
        else
            color = [.9 .2 .2 .25];
        end
        plot(Xs{i}(1:t,1),Xs{i}(1:t,2),'color',color)
    end
    xlim([0 2000])
    ylim([0 2000])
    %title(['t = ' num2str(round(T4(t),1))])
    ylabel('Migrants')
    xlabel('Non-migrants')

    set(gcf,'Position',[350 150 850 650])
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(findall(gcf,'-property','interpreter'),'interpreter','latex')

    %exportgraphics(gcf,'C:/Users/thekn/Pictures/Patritation/Phatot.gif',"Append",t~=1)
    %pause(.1)
end