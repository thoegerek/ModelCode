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

%%
[f,eqs,stability] = pipf(ax,ay,b,c,d,mu,5000,1/2.5);
%%
pip = struct;
 pip.data = {ax,ay,b,c,d,mu,X0};
 pip.dataDesc = {'ax','ay','b','c','d','mu','X0'};
 pip.tau = tau;
 pip.f = f;
 pip.eqs = eqs;
 pip.stability = stability;
 save('PIP.mat','pip');

% pip4.data = {ax,ay,b,c,d,mu,X0};
% pip4.dataDesc = {'ax','ay','b','c','d','mu','X0'};
% pip4.tau = tau;
% pip4.X = X;
% pip4.f = f;
% pip4.eqs = eqs;
% pip4.stability = stability;
% pips.d = .75;
% pips.ays = [.015;.0125;.015;.0125]
% pips.mus = [.24;.24;.22;.22]
% save('PIP123.mat','pip1','pip2','pip3','pip4','pips');
%%
figure(1)
imshow(((f(:,:,1)>0)*2+(f(:,:,2)>0))/3,'xdata',tau,'ydata',tau)
set(gca,'visible','on','ydir','normal')
xlabel('\tau_{r}')
ylabel('\tau_{i}')
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(eqs(stability<0),eqs(stability<0),'bo');
ueq = plot(eqs(stability>=0),eqs(stability>=0),'r-o');
legend([dummy1,dummy2,seq,ueq],{'(\partial \lambda*)/(\partial \tau) < 0','(\partial \lambda*)/(\partial \tau) > 0','stable equilibrium','unstable equilibrium'})
title('Pairwise invasibility plot')