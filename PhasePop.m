ax = .02;
ay = .01;
b = 1;
c = 1;
d = 1;
mu = .24;
tau = 0.1;
%%

[eqx,eqy,st] =equilibriumsStability(ax,ay,b,c,d,mu,tau,true);
steqx = eqx(st==-1);
steqy = eqy(st==-1);
usteqx = eqx(st~=-1);
usteqy = eqy(st~=-1);

Xspan = [0 2.5*max(steqx)];
Yspan = [0 2.5*max(steqy)+1];

%% Arrows
X = linspace(Xspan(1),Xspan(2),10);
Y = linspace(Yspan(1),Yspan(2),10);

nump1 = length(X);
nump2 = length(Y);

dx = zeros(nump1,nump2);
dy = zeros(nump1,nump2);

for i = 1:nump1
    for j = 1:nump2
        partial = myModel(0,[X(i);Y(j)],ax,ay,b,c,d,mu,tau);
        dx(j,i) = partial(1);
        dy(j,i) = partial(2);
    end
end
%% Transparrent lines
num2 = 25;%50;
rng(1)
X0r = rand(num2)*Xspan(2);
Y0r = rand(num2)*Yspan(2);
XT = cell(num2);
for i = 1:num2
    for j = 1:num2
        [~,XT{i,j}] = ode45(@myModel,Trange,[X0r(i,j);Y0r(i,j)],[],ax,ay,b,c,d,mu,tau);
    end
    disp([num2str(i) ' / ' num2str(num2) ' (2)'])
end
%% Solutions to define basins of attraction
Trange = (0:.02:10).^3;

num1 = 250;

X0 = [Xspan(1)*ones(1,num1),Xspan(2)*ones(1,num1)];
Y0 = repmat(linspace(Yspan(1),Yspan(2)+10,num1),1,2);
C = inf*ones(num1,1);
tol = 1e-2;
for i = 1:num1*2
    [~,Xt] = ode45(@myModel,Trange,[X0(i);Y0(i)],[],ax,ay,b,c,d,mu,tau);
    for j = 1:length(steqx)
        if (Xt(end,1)-steqx(j))^2 + (Xt(end,2)-steqy(j))^2 < tol
            C(i) = j;
        end
    end
    disp([num2str(i) ' / ' num2str(2*num1) ' (1)'])
end
%% Find the Seperatrix between the two basins
Saddle = [usteqx(usteqx~=0);usteqy(usteqx~=0)];
Sepind = find(diff(C)==1);
if length(Sepind) == 2
    SepY = (Y0(Sepind)+Y0(Sepind+1))/2;
    SepX = X0(Sepind);
    Sep = cell(2,1);
    for i = 1:2
        [~, Sep{i}] = ode45(@myModel,Trange,[SepX(i);SepY(i)],[],ax,ay,b,c,d,mu,tau);
        [~,SaddleMin] = min((Sep{i}(:,1)-Saddle(1)).^2+(Sep{i}(:,2)-Saddle(2)).^2);
        Sep{i} = Sep{i}(1:SaddleMin,:);
    end
    Seperatrix = [Sep{1};flipud(Sep{2})];
else
    Seperatrix = [nan nan];
end
%% Attracotrs
[~,A0] = ode15s(@myModel,Trange,[.0001,0],[],ax,ay,b,c,d,mu,tau);
[~,A1] = ode15s(@myModel,Trange,[10*Xspan(2),10*Yspan(2)],[],ax,ay,b,c,d,mu,tau);
if length(Saddle) > 1
    [~,A2] = ode15s(@myModel,Trange,[Saddle(1)-.001,Saddle(2)-.001],[],ax,ay,b,c,d,mu,tau);
    [~,A3] = ode15s(@myModel,Trange,[Saddle(1)+.0001,Saddle(2)+.0001],[],ax,ay,b,c,d,mu,tau);
else
    A2 = [nan nan];
    A3 = [nan nan];
end
%%
set(0,'defaulttextInterpreter','latex') 
figure(1)
hold on

if length(steqx)>1
    red = plot(polyshape([Seperatrix(:,1);Xspan(2);Xspan(1)],[Seperatrix(:,2);Yspan(1)-1;Yspan(1)-1]),'facecolor',[1 .4 .4],'edgealpha',0);
    green = plot(polyshape([Seperatrix(:,1);Xspan(2);Xspan(1)],[Seperatrix(:,2);Yspan(2)+1;Yspan(2)+1]),'facecolor',[.7 1 .7],'edgealpha',0);
    plot(Seperatrix(:,1),Seperatrix(:,2),'--k','linewidth',2)
else
    if steqy == 0
        red = plot(polyshape([Xspan(1) Xspan(1) Xspan(2) Xspan(2)],[Yspan(1) Yspan(2) Yspan(2) Yspan(1)]),'facecolor',[1 .4 .4],'edgealpha',0);
        green = [];
    else
        green = plot(polyshape([Xspan(1) Xspan(1) Xspan(2) Xspan(2)],[Yspan(1) Yspan(2) Yspan(2) Yspan(1)]),'facecolor',[.7 1 .7],'edgealpha',0);
        red = [];
    end
end

for i = 1:num2
    for j = 1:num2
        plot(XT{i,j}(:,1),XT{i,j}(:,2),'-','Color',[0,0,0,.01])
    end
end

sols = plot(nan,nan,'color',[.5 .5 .5]);
plot(A0(:,1),A0(:,2),'-k');
plot(A1(:,1),A1(:,2),'-k');
plot(A2(:,1),A2(:,2),'-k')
plot(A3(:,1),A3(:,2),'-k')

quiver(X,Y,dx,dy,'linewidth',1.5, 'AutoScaleFactor',0.5)    

psta = plot(steqx,steqy,'o','markersize',10,'linewidth',2,'color',[0.3 0.7 1]);
pusta = plot(usteqx,usteqy,'o','markersize',10,'linewidth',2,'color',[0.8500 0.3250 0.0980]);
xlabel('Non-migratory population ($n_{0}$)')
ylabel('Migratory population ($n_{1}$)')

legend([psta pusta,sols,red,green],{'Stable population equilibrium','Saddle point (Unstable eq. on manifold)','Solutions from random initial states','Converges to zero migrants','Converges to non-zero migrants'},'location','nw')
xlim([Xspan(1) Xspan(2)])
ylim([Yspan(1)-1 Yspan(2)+1])

set(gcf,'Position',[550 350 700 500])
%export_fig('C:/Users/thekn/Pictures/PhasePlotPopulationDynamics','-png','-transparent','-m5')