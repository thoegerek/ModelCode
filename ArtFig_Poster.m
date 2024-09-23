close all
%%
load('Invasion.mat')
load('PIP.mat')
load('PIP2.mat')
load('PIP3.mat')
load('PIP4.mat')
load('Bifurcation_tau.mat')
load('PhasePop.mat')
load('PiB_d.mat')
load('PiB_ay.mat')
load('PiB_mu.mat')
load('PopFromBif.mat')
load('PopFromBif_ay.mat')
load('PopFromBif_mu.mat')
load('Bistability_d.mat') 
load('Pop_Plot_1D.mat')
load('Pop_tau.mat')
% load('AXAY.mat')
% load('PiB_dimless.mat')
load('Pha.mat');
load('F_Bifurcation_tau.mat')
% load('Hys.mat');

% Color (and names) definitions
C = struct;
C.stable = [0.3 0.7 1];
C.unstable = [0.8500 0.3250 0.0980];
C.opt = [0.8500 0.3250 0.0980];%[.8 .2 .2];%[0.85 0.2 0.85];
C.totpop = [0 0 0];
C.nonmig = [0.6235    0.4000    0.7176];
C.mig = [0.25 0.65 0.1];
C.fraction = [0 0.4470 0.7410];
C.bireg = [.8 .2 .2];

names = cell(13,1);
%% Bifurcation in population eq.
figure
names{get(gcf,'number')} = 'Bifurcation_pop';

%stax1b
[stax,ord] = sort(bif.sx);
tstax = bif.tsx(ord);
stax1a = stax(tstax < .3);
tstax1a = tstax(tstax < .3);

stax0 = stax(logical((tstax >= .3 ).* (tstax <= .5)));
tstax0 = tstax(logical((tstax >= .3 ).* (tstax <= .5)));
%stax1b
tstax2a = tstax0(stax0 > 550);
stax2a = stax0(stax0 > 550);
[tstax2a,ord] = sort(tstax2a);
stax2a = stax2a(ord);

%stax2a
tstax1b = tstax0(stax0 <= 550);
stax1b = stax0(stax0 <= 550);

%stax2b
stax2b = stax(tstax > .5);
tstax2b = tstax(tstax > .5);

stax1 = [stax1a;stax1b];
tstax1 = [tstax1a,tstax1b];
stax2 = [stax2a;flip(stax2b)];
tstax2 = [tstax2a,flip(tstax2b)];

%ustax
[ustax,ord] = sort(bif.ux);
tustax = bif.tux(ord);
tustax = tustax(ustax>0);
ustax = ustax(ustax>0);

%stay1 
[stay,ord] = sort(bif.sy);

tstay = bif.tsy(ord);
tstay1 = tstay(stay==0);
stay1 = stay(stay==0);

%stay2 
tstay2 = tstay(stay~=0);
stay2 = stay(stay~=0);
[tstay2,ord] = sort(tstay2);
stay2 = stay2(ord);

%ustay
[ustay,ord] = sort(bif.uy);
tustay = bif.tuy(ord);
tustay = tustay(ustay>0);
ustay = ustay(ustay>0);
tustay = [max(tustax),tustay(end),tustay];
ustay = [0;0;ustay];

ymax = 1400;
ymin = 0; % maybe -1?


pat = patch([min(tustay) min(tustay) max(bif.tsy(bif.sy==0)) max(bif.tsy(bif.sy==0))],[ymin-1 ymax*1.1  ymax*1.1 ymin-1],C.bireg);
pat.EdgeColor = C.bireg;
pat.FaceAlpha = 0.2;
hold on

% saddlex1 = [min(tustax) (max(ustax)+min(stax2a))/2];
% saddley1 = [min(tustay) (max(ustay)+min(stay2))/2];
% saddlex2 = [max(tstax1)+(tau(2)-tau(1)) min(stax1)];
% saddley2 = [max(tstay1)+(tau(2)-tau(1)) 0];
% plot(saddlex1(1),saddlex1(2),'o',saddley1(1),saddley1(2),'o',saddlex2(1),saddlex2(2),'o',saddley2(1),saddley2(2),'o')

plt = plot(tustay,ustay,'--',tstay1,stay1,tstay2,stay2,tustax,ustax,'--',tstax1,stax1,tstax2,stax2,'linewidth',2);
colororder([
    C.mig;
    C.mig;
    C.mig;
    C.nonmig;
    C.nonmig;
    C.nonmig]);
dummy = plot(nan,nan,'k-',nan,nan,'k--','linewidth',2);
nothing = plot(nan,'color',[0 0 0 0]);

lgd = legend([plt(5),plt(2),dummy(1),dummy(2)],{'Non-migrants','Migrants','Stable equilibria','Unstable eq.'},'location','n');
lgd.NumColumns = 2;

xlim([bif.tau(1) bif.tau(end)])
ylim([ymin ymax])
xlabel('Sociality ($\sigma$)')
yticks([500 1000])
ylabel('Population')

set(gcf,'Position',[550 350 700 500])
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%%
axis
xlabel('time [years]')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
%%
path = 'C:/Users/thekn/Pictures/Poster/';

for i = 1:length(names)
    if ishandle(i)
        figure(i)
        export_fig([path,names{i}],'-png','-transparent','-m5')
    end
end