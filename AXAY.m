M = linearisedInvader(false);

Ax = 0.001:0.001:0.1;
Ay = 0.001:0.001:0.1;
b = 1;
c = .1;
d = .125;
mu = .5;

prec = 100;
uppercutoff = 1/3;
dtau = ((b-mu)/c)/(prec-1) *uppercutoff;
Tau = 0:dtau:(b-mu)/c  *uppercutoff;

%%
X = cell(length(Ax),length(Ay));
TotPop = zeros(length(Ax),length(Ay));
xPop = zeros(length(Ax),length(Ay));
yPop = zeros(length(Ax),length(Ay));

for i = 1:length(Ax)
    for j = 1:length(Ay)
    X0 = [(b-mu)/Ax(i),(b-mu)/Ay(j)];

    [eq,st] = evolutionaryEq(@myModel,X0,Ax(i),Ay(j),b,c,d,mu,Tau,0);
    tau = max(eq(st==-1));
    if isempty(tau)
        tau = 0;
    end
        [~,X{i,j}] = runToSS(@myModel,1,X0,1e3,1e-7,{Ax(i),Ay(j),b,c,d,mu,tau});

        TotPop(j,i) = sum(X{i,j}(end,:));
        xPop(j,i) = X{i,j}(end,1);
        yPop(j,i) = X{i,j}(end,2);

    disp([num2str((i-1)*length(Ax)+j) ' / ' num2str(length(Ax)*length(Ay))])
    end
end
prop = yPop./xPop;
prop(isnan(prop)) = 0;
%%
axay.data = {b,c,d,mu,X0};
axay.dataDesc = {'b','c','d','mu','X0'};
axay.ax = Ax;
axay.ay = Ay;
axay.TotPop = TotPop;
axay.prop = prop;
save('AXAY.mat','axay');