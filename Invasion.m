ax = .02;
ay = .015;
b = 1;
c = 1;
d = 1;
mu = .24;

taur = .2;
taui = .1;

%%
X01 = [1;1;0;0];
[T1,X1] = ode15s(@myModelInvader,[0 250],X01,[],ax,ay,b,c,d,mu,taur,0);

X02 = X1(end,:)' + [0;0;1;0];
[T2,X2] = ode15s(@myModelInvader,[0 500],X02,[],ax,ay,b,c,d,mu,taur,taui);

X03 = [X2(end,1); X2(end,2); X2(end,3); X2(end,4) - 15];
[T3,X3] = ode15s(@myModelInvader,[0 500],X03,[],ax,ay,b,c,d,mu,taur,taui);
%%
inv.data = {ax,ay,b,c,d,mu,taur,taui,X01,X02,X03};
inv.dataDesc = {'ax','ay','b','c','d','mu','taur','taui','X0','X02','X03'};
inv.T1 = T1;
inv.T2 = T2;
inv.T3 = T3;
inv.X1 = X1;
inv.X2 = X2;
inv.X3 = X3;
%save('Invasion.mat','inv');

%%
figure(1)
Tcat = [T1;T1(end)+T2;T1(end)+T2(end)+T3];
Xcat = [X1;X2;X3];
Xrcat = Xcat(:,1:2);
Xicat = Xcat(:,3:4);
ymax = 70;
plot([250 250],[0 ymax],'--',[750,750],[0,ymax],'--')
hold on
plt = plot(Tcat,Xrcat,Tcat,Xicat,'-.',Tcat,sum(Xcat,2),'linewidth',1.5);
colororder([.5 .5 .5;
    .5 .5 .5;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0 0.4470 0.7410]);
xlim([0,1000])
ylim([.5 ymax])
legend([plt(1),plt(2),plt(5)],{'Non-migrating population','Migrating population','Total population'},'position', [.35 .75 0 0] )
xlabel('Time (years)')
ylabel('Adult population number (1000s)')
