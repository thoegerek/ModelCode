a1 = .001;
a2 = .0003;
a3 = .0004;
a4 = .0005;
b = 1;
c = .2;
d = .015;
mu = .5;
tau = .4;

X0  = [420;100;0.1;0.001];
%%
[T,X] = ode45(@fourModel,[0 100],X0,[],a1,a2,a3,a4,b,c,d,mu,tau);
%%
figure(1)
plot(T,X)
legend('n1','n2','n3','n4')
%%
nump = 50;
N2 = linspace(0,100,nump);
N3 = linspace(0,100,nump);
Eq2 = zeros(nump);
Eq3 = zeros(nump);

for i = 1:nump
    for j = 1:nump
        X0 = [420;N2(i);N3(j);0];
        [~,X] = ode45(@fourModel,[0 1000],X0,[],a1,a2,a3,a4,b,c,d,mu,tau);
        Eq2(i,j) = sum(X(end,2));
        Eq3(i,j) = sum(X(end,3));
    end
    disp([num2str(i) ' / ' num2str(nump)])
end
%%
imagesc(Eq2+Eq3)
%% Multiple b_i
a1 = .001; 
a2 = .0003;
a3 = .0004;
a4 = .0005;
b1 = 1;
b2 = 1;
b2bad = .5;
b3 = 1;
b4 = 1;
c = .2;
d = .015;
mu = .5;
tau = .4;

X00 = [420;150;1;1];
%%
[T0,X0] = ode45(@fourModelMult,[0 100],X00,[],a1,a2,a3,a4,b1,b2,b3,b4,c,d,mu,tau);
X01 = X0(end,:);
[T1,X1] = ode45(@fourModelMult,[0 10],X01,[],a1,a2,a3,a4,b1,b2bad,b3,b4,c,d,mu,tau);
T1 = T1 + T0(end);
X02 = X1(end,:);
[T2,X2] = ode45(@fourModelMult,[0 100],X02,[],a1,a2,a3,a4,b1,b2,b3,b4,c,d,mu,tau);
T2 = T2 + T1(end);
%%
plot([T0;T1;T2],[X0;X1;X2])
hold on
plot([100 100],[0 max(X0,[],'all')*1.2],'--','color',[.7 .7 .7])
plot([110 110],[0 max(X0,[],'all')*1.2],'--','color',[.7 .7 .7])