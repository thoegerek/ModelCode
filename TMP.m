nump = 250;
mineps = .005;
maxeps = .1;
degree = 3;
Eps = (linspace(0,(maxeps-mineps)^(1/degree),nump).^degree) + mineps;

b = 1;
mu = .5;
del0 = .001;
del1 = .0003;
lambda = .2;

N0 = [100*(b-mu)/del0,100*(b-mu)];

prec1 = 150;
prec2 = 150;

tol = (1e-12)/(prec1+prec2)^2;
%% 1st grid
Eq = inf*ones(nump,1);
St = zeros(nump,1);

dsigma = ((b-mu)/lambda)/(prec1-1);
sigma = 0:dsigma:(b-mu)/lambda;

for i = 1:nump
    [eq,st] = evolutionaryEq(@myModel,N0,del0,del1,b,lambda,Eps(i),mu,sigma,tol);
    if length(eq) > size(Eq,2)
        Eq = [Eq,inf*ones(nump,1)];
        St = [St,inf*ones(nump,1)];
    end
    Eq(i,1:length(eq)) = eq;
    St(i,1:length(st)) = st;
    disp([num2str(i) ' / ' num2str(nump) ' (1)'])
end
%% 2nd grid
W = size(Eq,2);
WEq = inf*ones(nump*W,1);
WSt = 0*WEq;

for j = 1:W
    for i = 1:nump
        if Eq(i,j) < inf
        psigma = linspace(Eq(i,j)-dsigma*2,Eq(i,j)+dsigma*2,prec2);
        [eq,st] = evolutionaryEq(@myModel,N0,del0,del1,b,lambda,Eps(i),mu,psigma,tol);
        if length(eq) > size(WEq,2)
            WEq = [WEq,inf*ones(nump*W,1)];
            WSt = [WSt,inf*ones(nump*W,1)];
        end
        WEq(i+(j-1)*nump,1:length(eq)) = eq;
        WSt(i+(j-1)*nump,1:length(st)) = st;
        end
        disp([num2str(i+(j-1)*nump) ' / ' num2str(nump*W) ' (2)'])
    end
end

PEq = reshape(WEq,nump,[]);
PSt = reshape(WSt,nump,[]);
%% Data scrubbing that seems to be necessary due to numerical inaccuracies

PEq = reshape(WEq,nump,[]);
PSt = reshape(WSt,nump,[]);

bEq = true*ones(nump,size(PEq,2)-1);
for i = 1:nump
     [PEq(i,:),order] = sort(PEq(i,:));
     PSt(i,:) = PSt(i,order);
     bEq(i,:) = abs(diff(PEq(i,:))) > 5*dsigma/prec2;
     for j = 1:width(bEq)
         if ~bEq(i,j)
             PEq(i,j) = (PEq(i,j) + PEq(i,j+1))/2;
             PEq(i,j+1) = PEq(i,j);
         end
     end
end

CEq = PEq*0;
CSt = PSt*0;

for i = 1:nump
    [teq,order] = unique(PEq(i,:));
    tst = PSt(i,order);
    [teq,order] = sort(teq);
    CEq(i,:) = [teq, inf*ones(1,length(PEq(i,:))-length(teq))];
    CSt(i,:) = [tst(order),inf*ones(1,length(PEq(i,:))-length(teq))];
end
for i = size(CEq,2):-1:1
    if unique(CEq(:,i)) == inf
        CEq(:,i) = [];
        CSt(:,i) = [];
    end
end

%%
S_ESSeps = struct;
S_ESSeps.eq = CEq;
S_ESSeps.stability = CSt;
S_ESSeps.data = {del0,del1,b,lambda,mu,N0};
S_ESSeps.dataDesc = {'del0','del1','b','lambda','mu'};
S_ESSeps.eps = Eps;
save('FigD_Da.mat','S_ESSeps');
%%
figure(1)
axis;
hold on
for i = 1:size(CEq,2)
    plot(Eps(CSt(:,i)==1),CEq(CSt(:,i)==1,i),'m.')
    plot(Eps(CSt(:,i)==-1),CEq(CSt(:,i)==-1,i),'k.')
end
%%
Dsigma = .001;
Sigma = min(S_ESSeps.eq,[],'all'):Dsigma:max(S_ESSeps.eq(~isinf(S_ESSeps.eq)),[],'all')*10;
Bi = zeros(length(Sigma),length(S_ESSeps.d));
for i = 1:length(Eps)
    if sum(S_ESSeps.stability(i,:)==-1)>0
        onechance = false;
        taustart = find(Sigma-min(S_ESSeps.eq(i,:)) > 1/(prec1*prec2),1);
        for j = taustart:length(Sigma)
            [~,~,st] = equilibriumsStability(del0,del1,b,lambda,Eps(i),mu,Sigma(j),true);
            if sum(st==-1) > 1
                Bi(j,i) = 1;
            else
                if onechance
                    break;
                else
                    onechance = true;
                end
            end
            disp(['j = ' num2str(j) ' / ' num2str(length(Sigma)) ' (i = ' num2str(i) ' / ' num2str(length(S_ESSeps.d)) ' )'])
        end
    end
end
%%
pol = bwboundaries(Bi);
figure(2)
axis;
hold on
for i = 1:size(S_ESSeps.eq,2)
    plot(Eps(S_ESSeps.stability(:,i)==1),S_ESSeps.eq(S_ESSeps.stability(:,i)==1,i),'m.')
    plot(Eps(S_ESSeps.stability(:,i)==-1),S_ESSeps.eq(S_ESSeps.stability(:,i)==-1,i),'k.')
end
for i = 1:length(pol)
    patch(Eps(pol{i}(:,2)),Sigma(pol{i}(:,1)),'.r','facealpha',.2)
end
%%
S_eqeps.eq = S_ESSeps.eq;
S_eqeps.stability = S_ESSeps.stability;
S_eqeps.data = {S_ESSeps.data{1},S_ESSeps.data{2},S_ESSeps.data{3},S_ESSeps.data{4},S_ESSeps.data{5}};
S_eqeps.dataDesc = {'del0','del1','b','lambda','mu'};
S_eqeps.eps = S_ESSeps.eps;
S_eqeps.Sigma = Sigma;
S_eqeps.Bi = Bi;
save('FigD_Da.mat','S_eqeps',"-append");