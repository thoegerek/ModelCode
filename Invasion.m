ax = .001;
ay = .0003;
b = 1;
c = .2;
d = .015;
mu = .5;
mubad = .57;

tau1 = .8;
tau2 = .8;
tau3 = .4;
tau4 = .4;
%%
X00 = [0;0;1;1];
[T0,X0] = ode15s(@myModelInvader,[0 150],X00,[],ax,ay,b,c,d,mu,0,tau1);

X01 = [0;0;X0(end,3:4)'];
[T1,X1] = ode15s(@myModelInvader,[0 150],X01,[],ax,ay,b,c,d,mu,0,tau1);

X02 = [X1(end,1)+X1(end,3);X1(end,2)+X1(end,4);1;0];
[T2,X2] = ode15s(@myModelInvader,[0 600],X02,[],ax,ay,b,c,d,mu,tau1,tau2);
T2 = T2+T1(end);

X03 = [X2(end,1)+X2(end,3);X2(end,2)+X2(end,4);1;0];
[T3,X3] = ode15s(@myModelInvader,[0 600],X03,[],ax,ay,b,c,d,mu,tau2,tau3);
T3 = T3+T2(end);

X04 = [X3(end,1)+X3(end,3);X3(end,2)+X3(end,4);1;0];
[T4,X4] = ode15s(@myModelInvader,[0 1200],X04,[],ax,ay,b,c,d,mu,tau3,tau4);
T4 = T4+T3(end);


X05 = [X4(end,1)+X4(end,3);X4(end,2)+X4(end,4); 0 ; 0];
[T5,X5] = ode15s(@myModelInvader,[0 100],X05,[],ax,ay,b,c,d,mubad,tau4,0);
T5 = T5+T4(end);

X06 = [X5(end,1)+X5(end,3);X5(end,2)+X5(end,4); 0 ; 0];
[T6,X6] = ode15s(@myModelInvader,[0 150],X06,[],ax,ay,b,c,d,mu,tau4,0);
T6 = T6+T5(end);

Tcat = [T0-150;T1;T2;T3;T4;T5;T6];
Xcat = [X0;X1;X2;X3;X4;X5;X6];
%%
inv.data = {ax,ay,b,c,d,mu,mubad};
inv.dataDesc = {'ax','ay','b','c','d','mu','mubad'};
inv.taus = [tau1,tau2,tau3,tau4];
inv.T1 = T1;
inv.T2 = T2;
inv.T3 = T3;
inv.T4 = T4;
inv.T5 = T5;
inv.T6 = T6;
inv.X1 = X1;
inv.X2 = X2;
inv.X3 = X3;
inv.X4 = X4;
inv.X5 = X5;
inv.X6 = X6;
%save('Invasion.mat','inv');
%%
figure(1)
ymax = 2000;
plt = plot(T1,X1,T2,X2,T3,X3,T4,X4,T5,X5,T6,X6,Tcat,sum(Xcat,2),'linewidth',1.5);
colororder([...
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0 0.4470 0.7410]);
styles = {'-','-','-','-','-','-',...
'--','--','--','--',...
'-.','-.','-.','-.',...
':',':',...
':',':',':',':',':',':',':',':','-'};
lines = findobj(gca,'Type','line');
for i = 1:length(lines)
    lines(i).LineStyle = styles(length(styles)-i+1);
end
hold on
plot(T0-150,X0(:,4),'linewidth',1.5,'color',[0.9290 0.6940 0.1250])
plot(T0-150,X0(:,3),'linewidth',1.5,'color',[0.4660 0.6740 0.1880])
xlim([-150,Tcat(end)])
ylim([.5 ymax])
legend([plt(1),plt(2),plt(end)],{'Non-migrating population','Migrating population','Total population'},'location','ne')
xlabel('Time (years)')
ylabel('Adult population number (1000s)')
