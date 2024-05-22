load('PiB_d.mat')
Tau = bis.eq(:,2);
D = bis.d;
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
Xn = ((b-mu)/ax)*ones(length(D),1);
for i = 1:length(D)
    [eqx,eqy,stability] = equilibriumsStability(ax,ay,b,c,D(i),mu,Tau(i),false);
    X(i,:) = [max(eqx(stability==-1));max(eqy(stability==-1))];
    disp([num2str(i) ' / ' num2str(length(D))])
end
%%
pfb = struct;
pfb.X = X;
pfb.D = D;
pfb.data = {ax,ay,b,c,mu};
pfb.dataDesc = {'ax','ay','b','c','mu'};
save('PopFromBif.mat','pfb');
%%
figure(1)
plot(D,X(:,1),'g',D,X(:,2),'y')
hold on
plot(D,Xn(:,1),'g--')