D = [1:-.01:.05, 0.06:.01:1];
T = length(D);

ax = 0.01;
ay = 0.01;
b = 1;
c = .1;
%d = .1;
mu = .5;

uppercutoff = 1/4;
prec1 = 1000;
%%
X0 = [(b-mu)/ax,(b-mu)/ay];

x = zeros(T,2);
dtau = ((b-mu)/c)/(prec1-1) * uppercutoff;
tau = 0:dtau:(b-mu)/c *uppercutoff;
[eq,st] = evolutionaryEq(@myModel,X0,ax,ay,b,c,D(1),mu,tau,0);

Tau = zeros(T,1);

tauopt = max(eq(st==-1));


for t = 1:T
    [~,X] = runToSS(@myModel,1,X0,1e2,1e-10,{ax,ay,b,c,D(t),mu,tauopt});
    [eq,st] = evolutionaryEq(@myModel,X(end,:),ax,ay,b,c,D(t),mu,tau,0);
    if ~isempty(eq)
        diftau = eq-tauopt;
        pos = find(diftau >= 0);
        neg = find(diftau < 0);
        [over,id] = min(diftau(pos));
        ido = pos(id);
        [under,id] = max(diftau(neg));
        idu = neg(id);
        if ~isempty(over) && st(ido) == -1
            tauopt = over + tauopt;
        elseif ~isempty(under) && st(idu) == -1
            tauopt = under + tauopt;
        else
            tauopt = 0;
        end
    else
        tauopt = 0;
    end
    x(t,:) = X(end,:);
    Tau(t) = tauopt;

    X0 = x(t,:);

    disp([num2str(t) ' / ' num2str(T)])
end
%%
hys.tau = Tau;
hys.eq = x;
hys.data = {ax,ay,b,c,mu,X0};
hys.dataDesc = {'ax','ay','b','c','mu','X0'};
hys.D = D;
save('Hys.mat','hys');
%%
figure()
yyaxis left
plot(1:T,D,'k',1:T,Tau,'--b')
yyaxis right
plot(1:T,x)
legend('D','\tau*','x','y')