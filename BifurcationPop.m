% ax = .02;
% ay = .015;
% b = 1;
% c = 1;
% d = 1;
% mu = .24;

ax = .001;
ay = .0003;
b = 1;
c = .2;
d = .015;
mu = .5;


nump = 50;
cutoff = 1/3.5;
dtau = ((b-mu)/c)/(nump-1) * cutoff;
tau = 0:dtau:(b-mu)/c *cutoff;
%%
Eqx = zeros(length(tau),6);
Eqy = zeros(length(tau),6);
Q = zeros(length(tau),6);


for i = 1:length(tau)
    [eqx,eqy,stability] = equilibriumsStability(ax,ay,b,c,d,mu,tau(i),true);
    Eqx(i,1:length(eqx)) = eqx;
    Eqy(i,1:length(eqy)) = eqy;
    Q(i,1:length(stability)) = stability;
    disp([num2str(i) ' / ' num2str(length(tau))])
end
%%
taurep = repmat(tau,1,6);

sx = Eqx(Q==-1);
tsx = taurep(Q==-1);
ux = Eqx(Q~=-1);
tux = taurep(Q~=-1);

sy = Eqy(Q==-1);
tsy = taurep(Q==-1);
uy = Eqy(Q~=-1);
tuy = taurep(Q~=-1);
%%
bif.data = {ax,ay,b,c,d,mu};
bif.dataDesc = {'ax','ay','b','c','d','mu'};
bif.tau = tau;
bif.sx = sx;
bif.tsx = tsx;
bif.ux = ux;
bif.tux = tux;
bif.sy = sy;
bif.tsy = tsy;
bif.uy = uy;
bif.tuy = tuy;
%save('Bifurcation_tau.mat','bif');
%%
figure(3)
plot(tsx,sx,'.k',tux,ux,'.r',tuy,uy,'.m',tsy,sy,'.b');