load('PiB_d.mat')
Tau = pib.eq(:,2);
D = pib.d;
D(Tau == inf) = [];
Tau(Tau == inf) = [];
%%
ax = .001;
ay = .0003;
b = 1;
c = .2;
mu = .5;

%%
X = zeros(length(D),2);
xc = zeros(length(D),1);
U = zeros(length(D),2);
xn = ((b-mu)/ax)*ones(length(D),1);
for i = 1:length(D)
    [eqx,eqy,stability] = equilibriumsStability(ax,ay,b,c,D(i),mu,Tau(i),true);
    X(i,:) = [max(eqx(stability==-1));max(eqy(stability==-1))];
    if sum(stability==-1) > 1
        st = sort(eqx(stability==-1));
        xc(i) = st(1);
    end
    if sum(stability == 0) > 0
        if max(eqx(stability== 0)) > 0
            U(i,1) = max(eqx(stability== 0));
        end
        if max(eqy(stability ==  0)) > 0
            U(i,2) = max(eqy(stability ==  0));
        end
    end
    disp([num2str(i) ' / ' num2str(length(D))])
end
%%
pfb = struct;
pfb.X = X;
pfb.U = U;
pfb.xc = xc;
pfb.xn = xn;
pfb.D = D;


pfb.data = {ax,ay,b,c,mu};
pfb.dataDesc = {'ax','ay','b','c','mu'};
save('PopFromBif.mat','pfb');
%%
figure(1)
plot(D,X(:,1),'g',D,X(:,2),'r')
hold on
plot(D,xn(:,1),'g')
plot(D,U(:,1),'g--')
plot(D,U(:,2),'r--')
% plot(D,xc,'g')
%% for ay instead
load('PiB_ay.mat')
Tau = pibay.eq(:,2);
Ay = pibay.ay;
Ay(Tau == inf) = [];
Tau(Tau == inf) = [];
%%
ax = .001;
b = 1;
c = .2;
d = 0.015;
mu = .5;
%%
X = zeros(length(Ay),2);
xc = zeros(length(Ay),1);
U = zeros(length(Ay),2);
xn = ((b-mu)/ax)*ones(length(Ay),1);
for i = 1:length(Ay)
    [eqx,eqy,stability] = equilibriumsStability(ax,Ay(i),b,c,d,mu,Tau(i),true);
    X(i,:) = [max(eqx(stability==-1));max(eqy(stability==-1))];
    if sum(stability==-1) > 1
        st = sort(eqx(stability==-1));
        xc(i) = st(1);
    end
    if sum(stability == 0) > 0
        if max(eqx(stability== 0)) > 0
            U(i,1) = max(eqx(stability== 0));
        end
        if max(eqy(stability ==  0)) > 0
            U(i,2) = max(eqy(stability ==  0));
        end
    end
    disp([num2str(i) ' / ' num2str(length(Ay))])
end
%%
pfbay = struct;
pfbay.X = X;
pfbay.U = U;
pfbay.xc = xc;
pfbay.xn = xn;
pfbay.Ay = Ay;
pfbay.data = {ax,b,c,d,mu};
pfbay.dataDesc = {'ax','b','c','d','mu'};
%save('PopFromBif_ay.mat','pfbay');
%%
figure(1)
plot(Ay,X(:,1),'g')%,Ay,X(:,2),'r')
hold on
plot(Ay,xn(:,1),'g:')
plot(Ay,U(:,1),'g--')
plot(Ay,U(:,2),'r--')
plot(Ay,xc,'g.-')
%% For mu instead
load('PiB_mu.mat')
Tau = pibmu.eq(:,2);
Mu = pibmu.mu;
Mu(Tau == inf) = [];
Tau(Tau == inf) = [];
%%
ax = 10e-4;
ay = 3e-4;
b = 1;
c = .2;
d = 0.015;
%%
X = zeros(length(Mu),2);
xc = zeros(length(Mu),1);
U = zeros(length(Mu),2);
xn = ((b-Mu)/ax);
for i = 1:length(Mu)
    [eqx,eqy,stability] = equilibriumsStability(ax,ay,b,c,d,Mu(i),Tau(i),true);
    X(i,:) = [max(eqx(stability==-1));max(eqy(stability==-1))];
    if sum(stability==-1) > 1
        st = sort(eqx(stability==-1));
        xc(i) = st(1);
    end
    if sum(stability == 0) > 0
        if max(eqx(stability== 0)) > 0
            U(i,1) = max(eqx(stability== 0));
        end
        if max(eqy(stability ==  0)) > 0
            U(i,2) = max(eqy(stability ==  0));
        end
    end
    disp([num2str(i) ' / ' num2str(length(Mu))])
end
%%
pfbmu = struct;
pfbmu.X = X;
pfbmu.U = U;
pfbmu.xc = xc;
pfbmu.xn = xn;
pfbmu.Mu = Mu;
pfbmu.data = {ax,ay,b,c,d,mu};
pfbmu.dataDesc = {'ax','ay','b','c','d'};
save('PopFromBif_mu.mat','pfbmu');
%%
figure(1)
plot(Mu,X(:,1),'g',Mu,X(:,2),'r')
hold on
plot(Mu,U(:,1),'g--')

plot(Mu,xc,'g.-')
plot(Mu,xn,'k:')
plot(Mu,U(:,2),'r--')