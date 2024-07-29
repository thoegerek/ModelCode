M = linearisedInvaderComplete(true);

Ax = .001;
Ay = .0003;
B = 1;
C = .2;
D = .015;
Mu = .5;

%%
nump = 250;
f = zeros(nump,nump,2);
uppercutoff = 1/2.5;
dtau = ((B-Mu)/C)/(nump-1) * uppercutoff;
tau = 0:dtau:(B-Mu)/C *uppercutoff;
%%
X = zeros(nump,2);
Y = zeros(nump,2);
for r = 1:nump
    [eqx,eqy,stability] = equilibriumsStability(Ax,Ay,B,C,D,Mu,tau(r),true); %This is very slow - try to do this from (x0,0) instead of globally?
    nStab = sum(stability == -1);
    X(r,1:nStab) = eqx(stability == -1);
    Y(r,1:nStab) = eqy(stability == -1);

    [Y(r,1:nStab),ord] = sort(Y(r,1:nStab));
    X(r,1:nStab) = X(r,ord);
    
    for j = 1:2
        for i = 1:nump
            symM = subs(M,{'ax','ay','b','c','d','mu','taui','taur','xr','yr','xi','yi'},{Ax,Ay,B,C,D,Mu,tau(i),tau(r),X(r,j),Y(r,j),0,0});
            [~,lambda] = eig(double(symM));
            f(i,r,j) = max(diag(lambda));
        end
        if nStab == 1
            f(:,r,2) = f(:,r,1);
            break;
        end
    end

    disp([num2str(r) ' / ' num2str(nump)])
end
%%
figure(1)
imshow(((f(:,:,1)>0)*2+(f(:,:,2)>0))/3,'xdata',tau,'ydata',tau)
set(gca,'visible','on','ydir','normal')
xlabel('\tau_{r}')
ylabel('\tau_{i}')
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
legend([dummy1,dummy2],{'(\partial \lambda*)/(\partial \tau) < 0','(\partial \lambda*)/(\partial \tau) > 0'})
title('Pairwise invasibility plot')