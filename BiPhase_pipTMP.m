nump1 = 40;
minD = 0;
maxD = 5;
D = linspace(minD,maxD,nump1);

nump2 = 40;
minC = 1/nump2;
maxC = 2;
C = linspace(minC,maxC,nump2);

ax = 0.01;
ay = 0.01;
b = 1;
%c = .1;
%d = .1;
mu = .5;

X0 = [(b-mu)/ax,(b-mu)/ay];

maxtau = .1;

prec = 100;
tol = (1e-12)/prec^2;
%%
EQ = cell(nump2,1);
ST = cell(nump2,1);

for j = 1:nump2
    Eq = inf*ones(nump1,1);
    St = inf*ones(nump1,1);

    dtau = ((b-mu)/maxtau)/(prec-1);
    tau = 0:dtau:(b-mu)/maxtau;

    for i = 1:nump1
        [eq,st] = evolutionaryEq(@myModel,X0,ax,ay,b,C(j),D(i),mu,tau,tol,0);
        if length(eq) > size(Eq,2)
            Eq = [Eq,inf*ones(nump1,1)];
            St = [St,inf*ones(nump1,1)];
        end
        Eq(i,1:length(eq)) = eq;
        St(i,1:length(st)) = st;
        disp(['d bif. ' num2str(i) ' / ' num2str(nump1) ', c bif. ' num2str(j) ' / ' num2str(nump2)])
    end
    EQ{j} = Eq;
    ST{j} = St;
end
%%
B = zeros(nump1,nump2);

for i = 1:nump1
    for j = 1:nump2
        B(i,j) = sum((ST{j}(i,:)==-1));
    end
end
%%
figure
imshow(B==1,'xdata',C,'ydata',D)
set(gca,'visible','on','ydir','normal')
hold on
plot(C,C.^2*b/sqrt(ax))
xlabel('C')
ylabel('D')