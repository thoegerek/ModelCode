ax = .001;
ay = .0003;
b = 1;
c = .2;
d = .015;
mu = .5;
tau = .35; 
%%

[eqx,eqy,st] =equilibriumsStability(ax,ay,b,c,d,mu,tau,true);
steqx = eqx(st==-1);
steqy = eqy(st==-1);
usteqx = eqx(st~=-1);
usteqy = eqy(st~=-1);

Xspan = [0 3*max(steqx)];
Yspan = [0 3*max(steqy)];

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

%% Solutions to define basins of attraction
Trange = (0:.01:20).^2;

num1 = 1500;

X0 = [Xspan(1)*ones(1,num1),Xspan(2)*ones(1,num1)];
Y0 = repmat(linspace(Yspan(1),Yspan(2)+1,num1),1,2);
C = inf*ones(num1*2,1);
tol = 1e-0;
for i = 1:num1*2
    [~,Xt] = ode45(@myModel,Trange,[X0(i);Y0(i)],[],ax,ay,b,c,d,mu,tau);
    for j = 1:length(steqx)
        if (Xt(end,1)-steqx(j))^2 + (Xt(end,2)-steqy(j))^2 < tol
            C(i) = j;
        end
    end
    disp([num2str(i) ' / ' num2str(2*num1) ' (1)'])
end
%% Transparrent lines
num2 = 75;
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


%% Find the Seperatrix between the two basins
Saddle = [usteqx(usteqx~=0);usteqy(usteqy~=0)];
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
[~,A1] = ode15s(@myModel,Trange,[100,100],[],ax,ay,b,c,d,mu,tau);
if length(Sepind) == 2
    [~,A2] = ode15s(@myModel,Trange,[Saddle(1),Saddle(2)-.00001],[],ax,ay,b,c,d,mu,tau);
    [~,A3] = ode15s(@myModel,Trange,[Saddle(1),Saddle(2)+.000001],[],ax,ay,b,c,d,mu,tau);
else
    A2 = [nan nan];
    A3 = [nan nan];
end
%%
php = struct;
php.data = {ax,ay,b,c,d,mu,tau};
php.dataDesc = {'ax','ay','b','c','d','mu','tau'};
php.X = X;
php.Y = Y;
php.steqx = steqx;
php.steqy = steqy;
php.usteqx = usteqx;
php.usteqy = usteqy;
php.Xspan = Xspan;
php.Yspan = Yspan;
php.Seperatrix = Seperatrix;
php.XT = XT;
php.dx = dx;
php.dy = dy;
%save('PhasePop.mat','php');
%%
figure(1)

red = plot(polyshape([Seperatrix(:,1);Xspan(2);Xspan(1)],[Seperatrix(:,2);Yspan(1)-1;Yspan(1)-1]),'facecolor',[1 .7 .7],'edgealpha',0);
hold on
green = plot(polyshape([Seperatrix(:,1);Xspan(2);Xspan(1)],[Seperatrix(:,2);Yspan(2)+1;Yspan(2)+1]),'facecolor',[.7 1 .7],'edgealpha',0);
plot(Seperatrix(:,1),Seperatrix(:,2),'--k','linewidth',2)

for i = 1:num2
    for j = 1:num2
        plot(XT{i,j}(:,1),XT{i,j}(:,2),'-','Color',[0,0,0,.01])
    end
end

sols = plot(nan,nan,'color',[.5 .5 .5]);%plot(A1(:,1),A1(:,2),'-k');
%plot(A2(:,1),A2(:,2),'-k')
%plot(A3(:,1),A3(:,2),'-k')

quiver(X,Y,dx,dy,'linewidth',1.5, 'AutoScaleFactor',0.5)

psta = plot(steqx,steqy,'bo','markersize',10,'linewidth',2);
pusta = plot(usteqx,usteqy,'ro','markersize',10,'linewidth',2);
xlabel('Non-migratory population (n_{0})')
ylabel('Migratory population (n_{1})')
legend([psta pusta,sols,red,green],{'Stable population equilibrium','Saddle point (Unstable eq. on manifold)','Solutions from random initial states','Converges to zero mmigrators','Converges to non-zero migrators'},'location','nw')
xlim([Xspan(1) Xspan(2)])
ylim([Yspan(1)-.1 Yspan(2)+.1])

set(gcf,'Position',[550 350 700 500])