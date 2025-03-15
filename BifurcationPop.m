% ax = .001;
% ay = .0003;
% b = 1;
% c = .2;
% d = .015;
% mu = .5;

ax = .001;
ay = .0003;
b = 1;
c = .2;
d = .015;
mu = .5;
tau = .4;


% nump = 5000;
% cutoff = 1/2.5;
% dtau = ((b-mu)/c)/(nump-1) * cutoff;
% tau = 0:dtau:(b-mu)/c *cutoff;

nump = 500;
b = linspace(.5,1.5,nump);
%%
Eqx = zeros(length(tau),6);
Eqy = zeros(length(tau),6);
Q = zeros(length(tau),6);


for i = 1:length(b)
    [eqx,eqy,stability] = equilibriumsStability(ax,ay,b(i),c,d,mu,tau,true);
    Eqx(i,1:length(eqx)) = eqx;
    Eqy(i,1:length(eqy)) = eqy;
    Q(i,1:length(stability)) = stability;
    disp([num2str(i) ' / ' num2str(length(b))])
end
%%
brep = repmat(b,1,6);

sx = Eqx(Q==-1);
tsx = brep(Q==-1);
ux = Eqx(Q~=-1);
tux = brep(Q~=-1);

sy = Eqy(Q==-1);
tsy = brep(Q==-1);
uy = Eqy(Q~=-1);
tuy = brep(Q~=-1);
%%
bif.data = {ax,ay,b,c,d,tau};
bif.dataDesc = {'ax','ay','b','c','d','tau'};
bif.mu = mu;
bif.sx = sx;
bif.tsx = tsx;
bif.ux = ux;
bif.tux = tux;
bif.sy = sy;
bif.tsy = tsy;
bif.uy = uy;
bif.tuy = tuy;
%save('Bifurcation_mu.mat','bif');
%%
figure(3)
plot(tsx,sx,'.k',tux,ux,'.r',tuy,uy,'.m',tsy,sy,'.b');