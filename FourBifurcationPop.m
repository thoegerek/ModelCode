a1 = .001;
a2 = .0006;
a3 = .0009;
a4 = .0018;
b = 1;
c = .2;
d = .015;
mu = .5;


nump = 2500;
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
syms A B C D MU TAU
LHS = MU/(D*TAU^2);
RHS = ((B-MU-C*TAU)^2)/A;

ss = solve(LHS==RHS,TAU);
UpperBound = double(subs(ss(4),{A,B,C,D,MU},{a1,b,c,d,mu})); 
clear A B C D MU TAU
%%
max1 = max(Eq1,[],2);
max2 = max(Eq2,[],2);
max3 = max(Eq3,[],2);
max4 = max(Eq4,[],2);
csum = max1.^2*a1 + max2.^2*a2 + max3.^2*a3 + max4.^2*a4;
mls = mu./(d*tau.^2);
[~,lower] = min(abs(csum-mls'));
LowerBound = tau(lower);
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
f_bif.ub = UpperBound;
f_bif.lb = LowerBound;
f_bif.sr = max2+max3+max4;
f_bif.mls = mls;
f_bif.csum = csum;
f_bif.thresh = (b-mu-c*tau).^2/a1;
%save('F_Bifurcation_tau.mat','f_bif');

%%
figure(3)
ymax = max([Eq1 Eq2 Eq3 Eq4],[],'all')*1.2;
plot([UpperBound,UpperBound],[0 ymax],'--k',[LowerBound,LowerBound],[0 ymax],'--k')
hold on
plot(ts,s1,'*k',tu,u1,'.k',tu,u2,'.b',ts,s2,'*b',tu,u3,'.g',ts,s3,'*g',tu,u4,'.r',ts,s4,'*r','markersize',10);
legend('n1','n2','n3','n4')
ylim([0 ymax])