nump1 = 5;
minMu = 1/(nump1+1);
maxMu = 1-1/nump1;
Mu = linspace(minMu,maxMu,nump1);

nump2 = 5;
minlogAy = -6;
maxlogAy = 6;
logAy = linspace(minlogAy,maxlogAy,nump2);
Ay = exp(logAy);


D = kron(10.^(-5:5),[1,3]);
nump3 = length(D);

prec = 1000;
tol = (1e-8)/prec^2;
%%
tic
A = zeros(nump1,nump2);
B = zeros(nump1,nump2);
C = zeros(nump1,nump2);
for i = 1:nump1
    for j = 1:nump2
        currprec = prec;
        X0 = [1-Mu(i),(1-Mu(i))/Ay(j)];
        degree = 1;%+1.5^logAy(j); %Houristic for hitting all tau posibilities while still keeping prec down
        tau = linspace(0,(1-Mu(i))^(1/degree),currprec).^degree;
        k = max(length(D),B(i,max(1,j-1)));
        while k >= 1 && B(i,j) == 0
            [eq,st] = evolutionaryEq(@myModel,X0,1,Ay(j),1,1,D(k),Mu(i),tau,tol);
            if ~isempty(eq(st==-1)>0)
                A(i,j) = D(k);
                B(i,j) = k;
                C(i,j) = max(eq);
                k = k-1;
            end
            if B(i,j) == 0
                break;
%                 currprec = currprec*3;
%                 tau = linspace(0,(1-Mu(i)),currprec);
%                 k = max(length(D),B(i,max(1,j-1)));
%                 disp(['Failed to find value, increasing precision to ' num2str(currprec)])
            else
                disp(['ay: ' num2str(j) ' / ' num2str(nump2) ', mu: ' num2str(i) ' / ' num2str(nump1)])
            end
        end
    end
end
toc
%%
% E = zeros(nump1,nump2);
% prec2 = 1000;
% for i = 1:nump1
%     for j = 1:nump2
%         X0 = [1-Mu(i),(1-Mu(i))/Ay(j)];
%         tau = linspace(0,C(i,j),prec2);
%         for k = B(i,j):length(D)
%             [eq,st] = evolutionaryEq(@myModel,X0,1,Ay(j),1,1,D(k),Mu(i),tau,tol);
%             if ~isempty(eq(st==-1)>0)
%                 [eqx,~,stx] = equilibriumsStability(1,Ay(j),1,1,D(k),Mu(i),eq(st==-1),true);
%                 if length(eqx(stx==-1)) > 1
%                     E(i,j) = eq(st==-1);
%                 end
%             end
%         end
%         disp(['ay: ' num2str(j) ' / ' num2str(nump2) ', mu: ' num2str(i) ' / ' num2str(nump1)])
%     end
% end
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