a1 = .001;
a2 = .0003;
a3 = .0004;
a4 = .0005;
b = 1;
c = .2;
d = .015;
mu = .5;


nump = 5000;
cutoff = 1/2.5;
dtau = ((b-mu)/c)/(nump-1) * cutoff;
tau = 0:dtau:(b-mu)/c *cutoff;
%%
Eq1 = zeros(length(tau),30);
Eq2 = zeros(length(tau),30);
Eq3 = zeros(length(tau),30);
Eq4 = zeros(length(tau),30);
Q = zeros(length(tau),30);


for i = 1:length(tau)
    [eq1,eq2,eq3,eq4,stability] = fourEquilibriumsStability(a1,a2,a3,a4,b,c,d,mu,tau(i),true);
    Eq1(i,1:length(eq1)) = eq1;
    Eq2(i,1:length(eq2)) = eq2;
    Eq3(i,1:length(eq3)) = eq3;
    Eq4(i,1:length(eq4)) = eq4;
    Q(i,1:length(stability)) = stability;
    disp([num2str(i) ' / ' num2str(length(tau))])
end
%%
taurep = repmat(tau,1,30);

s1 = Eq1(Q==-1);
u1 = Eq1(Q~=-1);

s2 = Eq2(Q==-1);
u2 = Eq2(Q~=-1);

s3 = Eq3(Q==-1);
u3 = Eq3(Q~=-1);

s4 = Eq4(Q==-1);
u4 = Eq4(Q~=-1);

ts = taurep(Q==-1);
tu = taurep(Q~=-1);
%%
f_big = struct;
f_bif.data = {a1,a2,a3,a4,b,c,d,mu};
f_bif.dataDesc = {'a1','a2','a3','a4','b','c','d','mu'};
f_bif.tau = tau;
f_bif.s1 = s1;
f_bif.u1 = u1;
f_bif.s2 = s2;
f_bif.u2 = u2;
f_bif.s3 = s3;
f_bif.u3 = u3;
f_bif.s4 = s4;
f_bif.u4 = u4;
f_bif.ts = ts;
f_bif.tu = tu;
%save('F_Bifurcation_tau.mat','f_bif');
%%
figure(3)
plot(ts,s1,'*k',tu,u1,'.k',tu,u2,'.b',ts,s2,'*b',tu,u3,'.g',ts,s3,'*g',tu,u4,'.r',ts,s4,'*r','markersize',1);
legend('n1','n2','n3','n4')