nump = 250;
minAy = 0.00002;
maxAy = 0.0007;
Ay = linspace(minAy,maxAy,nump);

ax = .001;
ay = .0003;
b = 1;
c = .2;
d = .015;
mu = .5;

X0 = [100*(b-mu)/ax,100*(b-mu)];

uppercutoff = 1;
prec1 = 150;
prec2 = 150;

tol = (1e-12)/(prec1+prec2)^2;
%% 1st grid
Eq = inf*ones(nump,1);
St = zeros(nump,1);

dtau = ((b-mu)/c)/(prec1-1) * uppercutoff;
tau = 0:dtau:(b-mu)/c *uppercutoff;

for i = 1:nump
    [eq,st] = evolutionaryEq(@myModel,X0,ax,Ay(i),b,c,d,mu,tau,tol);
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
        ptau = linspace(Eq(i,j)-dtau*2,Eq(i,j)+dtau*2,prec2);
        [eq,st] = evolutionaryEq(@myModel,X0,ax,Ay(i),b,c,d,mu,ptau,tol);
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
     bEq(i,:) = abs(diff(PEq(i,:))) > 5*dtau/prec2;
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
pibay.eq = CEq;
pibay.stability = CSt;
pibay.data = {ax,ay,b,c,mu,X0};
pibay.dataDesc = {'ax','ay','b','c','mu','X0'};
pibay.ay = Ay;
%save('PiB_ay.mat','pibay');
%%
figure(1)
axis;
hold on
for i = 1:size(CEq,2)
    plot(Ay(CSt(:,i)==1),CEq(CSt(:,i)==1,i),'m.')
    plot(Ay(CSt(:,i)==-1),CEq(CSt(:,i)==-1,i),'k.')
end
%%
% load('PiB_ay.mat')
% %%
% Dtau = .001;
% Tau = min(pibay.eq,[],'all'):Dtau:max(pibay.eq(~isinf(pibay.eq)),[],'all')*10;
% Bi = zeros(length(Tau),length(pibay.ay));
% for i = 1:length(Ay)
%     if sum(pibay.stability(i,:)==-1)>0
%         onechance = false;
%         taustart = find(Tau-min(pibay.eq(i,:)) > 1/(prec1*prec2),1);
%         for j = taustart:length(Tau)
%             [~,~,st] = equilibriumsStability(ax,Ay(i),b,c,d,mu,Tau(j),true);
%             if sum(st==-1) > 1
%                 Bi(j,i) = 1;
%             else
%                 if onechance
%                     break;
%                 else
%                     onechance = true;
%                 end
%             end
%             disp(['j = ' num2str(j) ' / ' num2str(length(Tau)) ' (i = ' num2str(i) ' / ' num2str(length(pibay.ay)) ' )'])
%         end
%     end
% end
% %%
% bisay.eq = pibay.eq;
% bisay.stability = pibay.stability;
% bisay.data = {pibay.data{1},pibay.data{2},pibay.data{3},pibay.data{4},pibay.data{5}};
% bisay.dataDesc = {'ax','ay','b','c','mu'};
% bisay.ay = pibay.ay;
% bisay.Tau = Tau;
% bisay.Bi = Bi;
% save('Bistability_ay.mat','bisay');
% %%
% pol = bwboundaries(Bi);
% figure(2)
% axis;
% hold on
% for i = 1:size(pibay.eq,2)
%     plot(Ay(pibay.stability(:,i)==1),pibay.eq(pibay.stability(:,i)==1,i),'m.')
%     plot(Ay(pibay.stability(:,i)==-1),pibay.eq(pibay.stability(:,i)==-1,i),'k.')
% end
% for i = 1:length(pol)
%     patch(Ay(pol{i}(:,2)),Tau(pol{i}(:,1)),'.r','facealpha',.2)
% end