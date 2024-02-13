nump1 = 5;
minMu = 1/(nump1+1);
maxMu = 1-1/nump1;
Mu = linspace(minMu,maxMu,nump1);

nump2 = 5;
minlogAy = -6;
maxlogAy = 0;%6;
logAy = linspace(minlogAy,maxlogAy,nump2);
Ay = exp(logAy);


D = kron([1,3],10.^(-5:5));

prec = 300;
tol = (1e-16)/prec^2;
%%
tic
A = zeros(nump1,nump2);%ciritcal value
B = length(D)*ones(nump1,nump2);%index of ciritcal value w.r.t. D
C = zeros(nump1,nump2);
E = zeros(nump1,nump2);
for i = 1:nump1
    for j = 1:nump2
        X0 = [1-Mu(i),(1-Mu(i))/Ay(j)];
        tau = linspace(0,sqrt(1-Mu(i)),prec);%.^2; %define grid (denser near 0) of possible tau values
        prevst = [];
        for k = B(i,min(nump2,j+1)):-1:1 %max(1,B(max(1,i-1),j)-2) goes from critical value of the pixel below (-2 for safety), until the max value (or the break)
            [eq,st] = evolutionaryEq(@myModel,X0,1,Ay(j),1,1,D(k),Mu(i),tau,tol); %find eqs. for k
            if ~isempty(eq(st==-1)>0) %If a stable interioir eq. exists, write k as the critical value
                A(i,j) = D(k);
                B(i,j) = k;
                prevst = st;
            else
                steq = prevst==-1;
                usteq = prevst~=-1;
                if ~isempty(steq)
                    [~,eqy,~] = equilibriumsStability(1,Ay(j),1,1,D(max(k-1,1)),Mu(i),steq(1),true);
                    if length(eqy) > 2
                        C(i,j) = 1;
                    end
                end
                if ~isempty(usteq)
                    [~,eqy,~] = equilibriumsStability(1,Ay(j),1,1,D(max(k-1,1)),Mu(i),usteq(1),true);
                    if length(eqy) > 2
                        E(i,j) = 1;
                    end
                end
                break;
            end
        end
        disp(['ay: ' num2str(j) ' / ' num2str(nump2) ', mu: ' num2str(i) ' / ' num2str(nump1)])
    end
end
toc
%%
pha.Mu = Mu;
pha.Ay = Ay;
pha.D = D;
pha.bif = A;
pha.Dind = B;
pha.numparms = {prec,tol};
%save('Pha.mat','pha');
%%
figure()
imagesc(log(A),'xdata',logAy,'ydata',Mu)
xlabel('ln(n_{1}/n_{0})')
ylabel('\mu/b')
set(gca,'ydir','normal')