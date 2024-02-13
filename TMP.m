ax = 0.01;
ay = 0.01;
b = 1;
mu = .5;

nump1 = 10;
minD = .05;
maxD = 1;
D = linspace(minD,maxD,nump1);

nump2 = 10;
maxC = .1;
minC = maxC/nump2;
C = linspace(minC,maxC,nump2);

prec = 50;
maxtau = 1;
%%
ny = cell(nump1,nump2);
my = cell(nump1,nump2);

for i = 1:nump1
    for j = 1:nump2
        Tau = linspace(0,.25*(b-mu)/maxtau,prec);
        ny{i,j} = zeros(prec,1);
        for k = 1:prec
            [eqx,eqy,st] = equilibriumsStability(ax,ay,b,C(j),D(i),mu,Tau(k),false);
            ny{i,j}(k) = sum((st==-1));
            my{i,j}(k) = sum((eqy>0).*(st==-1));

            if mod(k,5) == 0
                disp(['tau: ' num2str(k) ' / ' num2str(prec) ', C: ' num2str(j) ' / ' num2str(nump2) ', D: ' num2str(i) ' / ' num2str(nump1)])
            end
        end
    end
end
%%
A = zeros(nump1,nump2);
for i = 1:nump1
    for j = 1:nump2
        A(i,j) = sum(ny{i,j});
        B(i,j) = sum(my{i,j});
    end
end