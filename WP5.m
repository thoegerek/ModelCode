%%
plot(Tcat,Xcat)
colororder([...
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.9290 0.6940 0.1250]);
%styles = {'-','-','-','-','-','-',...
%'--','--','--','--',...
%'-.','-.','-.','-.',...
%':',':',...
%':',':',':',':',':',':',':',':','-'};
%lines = findobj(gca,'Type','line');

%%
ax = .001;
ay = .0003;
b = 1;
c = .2;
d = .015;
mu = .5;
bbad = .9;
bworse = .7;

tau = .6;
%%
X00 = [(b-mu-c*tau)/ax;(b-mu-c*tau)/ay];
[T0,X0] = ode15s(@myModel,[0 100],X00,[],ax,ay,b,c,d,mu,tau);
T0 = T0-T0(end);

X10 = X0(end,:)'+[0;.01];
[T1,X1] = ode15s(@myModel,[0 50],X10,[],ax,ay,b,c,d,mu,tau);

X20 = X1(end,:)'+[0;.01];
[T2,X2] = ode15s(@myModel,[0 100],X20,[],ax,ay,bbad,c,d,mu,tau);
T2 = T2+T1(end);

X30 = X2(end,:)'+[0;.01];
[T3,X3] = ode15s(@myModel,[0 75],X30,[],ax,ay,bworse,c,d,mu,tau);
T3 = T3+T2(end);

X40 = X3(end,:)'+[0;.01];
[T4,X4] = ode15s(@myModel,[0 150],X40,[],ax,ay,bbad,c,d,mu,tau);
T4 = T4+T3(end);

X50 = X4(end,:)'+[0;.01];
[T5,X5] = ode15s(@myModel,[0 200],X50,[],ax,ay,b,c,d,mu,tau);
T5 = T5+T4(end);

Tcat = [T0;T1;T2;T3;T4;T5];
Xcat = [X0;X1;X2;X3;X4;X5];
%%
ymax = max(Xcat,[],'all')*1.2;
hold on
patch([T2(1) T2(end) T2(end) T2(1)], [ymax ymax 0 0]+5, 'r','facecolor', [.8 .8 .8], 'edgealpha', 0)
patch([T3(1) T3(end) T3(end) T3(1)], [ymax ymax 0 0]+5, 'r','facecolor', [.6 .6 .6], 'edgealpha', 0)
patch([T4(1) T4(end) T4(end) T4(1)], [ymax ymax 0 0]+5, 'r','facecolor', [.8 .8 .8], 'edgealpha', 0)
plot(Tcat,Xcat)
xlim([0,max(Tcat)])
ylim([0 ymax])
%%
