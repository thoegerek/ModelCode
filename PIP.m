%Main (small d small ay)
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

cutoff = 1/((b-mu)/c);
tau = linspace(0,cutoff*(b-mu)/c,nump);
%%
nump = 50;
[f,eqs,stability] = pipf(ax,ay,b,c,d,mu,nump,cutoff);
%%
pip = struct;
 pip.data = {ax,ay,b,c,d,mu};
 pip.dataDesc = {'ax','ay','b','c','d','mu'};
 pip.tau = tau;
 pip.f = f;
 pip.eqs = eqs;
 pip.stability = stability;
%save('PIP.mat','pip');
%%
figure(1)
imshow(((pip.f(:,:,1)>0)*2+(pip.f(:,:,2)>0))/3,'xdata',pip.tau,'ydata',pip.tau)
set(gca,'visible','on','ydir','normal')
xlabel('\tau_{r}')
ylabel('\tau_{i}')
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(pip.eqs(pip.stability<0),pip.eqs(pip.stability<0),'bo');
ueq = plot(pip.eqs(pip.stability>=0),pip.eqs(pip.stability>=0),'r-o');
%% large d large ay
ax = .001;
ay = .00075;
b = 1;
c = .2;
d = .025;
mu = .5;

cutoff = 1/((b-mu)/c);
tau = linspace(0,cutoff*(b-mu)/c,nump);
%%
nump = 2500;
[f,eqs,stability] = pipf(ax,ay,b,c,d,mu,nump,cutoff);
%%
pip2 = struct;
 pip2.data = {ax,ay,b,c,d,mu};
 pip2.dataDesc = {'ax','ay','b','c','d','mu'};
 pip2.tau = tau;
 pip2.f = f;
 pip2.eqs = eqs;
 pip2.stability = stability;
save('PIP2.mat','pip2');

%%
figure(2)
imshow(((pip2.f(:,:,1)>0)*2+(pip2.f(:,:,2)>0))/3,'xdata',pip2.tau,'ydata',pip2.tau)
set(gca,'visible','on','ydir','normal')
xlabel('\tau_{r}')
ylabel('\tau_{i}')
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(pip2.eqs(pip2.stability<0),pip2.eqs(pip2.stability<0),'bo');
ueq = plot(pip2.eqs(pip2.stability>=0),pip2.eqs(pip2.stability>=0),'r-o');
%% Large d and small ay
ax = .001;
ay = .0003;
b = 1;
c = .05;
d = .025;
mu = .5;

cutoff = 1/((b-mu)/c);
tau = linspace(0,cutoff*(b-mu)/c,nump);
%%
nump = 2500;
[f,eqs,stability] = pipf(ax,ay,b,c,d,mu,nump,cutoff);
%%
pip3 = struct;
 pip3.data = {ax,ay,b,c,d,mu};
 pip3.dataDesc = {'ax','ay','b','c','d','mu'};
 pip3.tau = tau;
 pip3.f = f;
 pip3.eqs = eqs;
 pip3.stability = stability;
save('PIP3.mat','pip3');
%%
figure(3)
imshow(((pip3.f(:,:,1)>0)*2+(pip3.f(:,:,2)>0))/3,'xdata',pip3.tau,'ydata',pip3.tau)
set(gca,'visible','on','ydir','normal')
xlabel('\tau_{r}')
ylabel('\tau_{i}')
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(pip3.eqs(pip3.stability<0),pip3.eqs(pip3.stability<0),'bo');
ueq = plot(pip3.eqs(pip3.stability>=0),pip3.eqs(pip3.stability>=0),'r-o');
%% Small d, Large ay
ax = .001;
ay = .00075;
b = 1;
c = .2;
d = .015;
mu = .5;

cutoff = 1/((b-mu)/c);
tau = linspace(0,cutoff*(b-mu)/c,nump);
%%
nump = 2500;
[f,eqs,stability] = pipf(ax,ay,b,c,d,mu,nump,cutoff);
%%
pip4 = struct;
 pip4.data = {ax,ay,b,c,d,mu};
 pip4.dataDesc = {'ax','ay','b','c','d','mu'};
 pip4.tau = tau;
 pip4.f = f;
 pip4.eqs = eqs;
 pip4.stability = stability;
save('PIP4.mat','pip4');
%%
figure(4)
imshow(((pip4.f(:,:,1)>0)*2+(pip4.f(:,:,2)>0))/3,'xdata',pip4.tau,'ydata',pip4.tau)
set(gca,'visible','on','ydir','normal')
xlabel('\tau_{r}')
ylabel('\tau_{i}')
hold on
dummy1 = plot(inf,inf,'squarek','markerfacecolor','k');
dummy2 = plot(inf,inf,'squarek');
seq = plot(pip4.eqs(pip4.stability<0),pip4.eqs(pip4.stability<0),'bo');
ueq = plot(pip4.eqs(pip4.stability>=0),pip4.eqs(pip4.stability>=0),'r-o');