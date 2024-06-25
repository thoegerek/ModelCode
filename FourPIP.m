a1 = .001;
a2 = .0003;
a3 = .0004;
a4 = .0005;
b = 1;
c = .2;
d = .015;
mu = .5;
%%
nump = 150;

cutoff = 1/((b-mu)/c);
tau = linspace(0,cutoff*(b-mu)/c,nump);
[f,eqs,stability] = fourPipf(a1,a2,a3,a4,b,c,d,mu,nump,cutoff);
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