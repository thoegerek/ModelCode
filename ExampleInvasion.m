ax = 0.02; %0.01
ay = 0.002;
b = 20;
c = 2;
d = .00125;
mu = 10;

taur = 1.25;
taui = .95;

%% No invader
X0 = [300;300];
[Tl,Xl] = ode45(@myModel,[0 5],X0,[],ax,ay,b,c,d,mu,taur);
X0 = [250;250];
[Td,Xd] = ode45(@myModel,[0 5],X0,[],ax,ay,b,c,d,mu,taur);
%%
figure(1)
plt = plot(Tl,Xl(:,1),Tl,Xl(:,2),'linewidth',1.5);
set(gca,'colororder',[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[0.8500 0.3250 0.0980]])
hold on
plot(Td,Xd(:,1),'--',Td,Xd(:,2),'--','linewidth',1.5)
dummy = plot([0 0],[-1 -1],'-',[0 0],[-1 -1],'--','color',[.1 .1 .1]);
legend([plt(1),plt(2),dummy(1),dummy(2)],{'x','y','surviving migratory route','collapsing migratory route'},'location','nw')
ylim([0 max(Xl,[],'all')*1.25])
ylabel('Adult individuals')
xlabel('Time [years]')
%% With invader
[eqx,eqy,st] = equilibriumsStability(ax,ay,b,c,d,mu,taur,true)

X0 = [0;0;1;0];
X0(1) = eqx(6);
X0(2) = eqy(6); % "6" found manually
%%
[T,X] = ode15s(@myModelInvader,[0 100],X0,[],ax,ay,b,c,d,mu,taur,taui);
%%
inv.data = {ax,ay,b,c,d,mu,X0};
inv.dataDesc = {'ax','ay','b','c','d','mu','X0'};
inv.taur = taur;
inv.taui = taui;
inv.T = T;
inv.X = X;
%save('Invasion.mat','inv');
%%
figure(1)
plot(T,X(:,1),T,X(:,2),T,X(:,3),':',T,X(:,4),':','linewidth',1.5);
set(gca,'colororder',[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[0.8500 0.3250 0.0980]])
hold on
plot(T,sum(X,2),'k--','linewidth',1.5)
legend('x_{\tau = 1.25}','y_{\tau = 1.25}','x_{\tau = .95}','y_{\tau = .95}','Total population')
ylabel('Adult individuals')
xlabel('Time [years]')
