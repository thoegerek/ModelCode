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