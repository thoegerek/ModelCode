nump1 = 5;
minlogB = .5;
maxlogB = 2;
logB = linspace(minlogB,maxlogB,nump1);
B = exp(logB);

nump2 = 6;
minlogAy = -6;
maxlogAy = 6;
logAy = linspace(minlogAy,maxlogAy,nump2);
Ay = exp(logAy);


D = kron(10.^(-5:10),[1,2,3,5]);
Dnan = [nan,D]; %For disping
nump3 = length(D);


%%
prec = [1, 10, 100, 1000, 10000];
A = zeros(nump1,nump2);%Index for D as critical value
P = zeros(nump1,nump2);%Precision utilised
TAU = zeros(nump1,nump2);%ESS
for i = 1:nump1 %loops over b
    for j = 1:nump2 %loops over Ay
        X0 = [B(i)-1,(B(i)-1)/Ay(j)]; %Apparently it's importnt this does not change
        for p = 1:length(prec)-1 %Loop that increases precision (and zooms) when no ESS is found
            tau = linspace(0,(B(i)-1)/prec(p),(prec(p+1))/prec(p));
            for k = 1:nump3
                [eq,st] = evolutionaryEq(@myModel,X0,1,Ay(j),B(i),1,D(k),1,tau,1e-8);
                if ~isempty(eq)
                    A(i,j) = k;
                    break  
                end
            end
            if A(i,j) > 0 %If A>0, then a value was found
                P(i,j) = prec(p);
                TAU(i,j) = eq(end);
                if j == 1 || TAU(i,j) <= TAU(i,j-1) %Checks for false positives arising when D large and Ay large [REFINE THIS!]
                    break
                else %If so, throw away answer and increase precision
                    A(i,j) = 0;
                    P(i,j) = 0;
                    TAU(i,j) = 0;
                end
            end
        end
        disp(['D = ' num2str(Dnan(A(i,j)+1)) ' (' num2str(i) '/' num2str(nump1) ') (' num2str(j) '/' num2str(nump2) ')'])
    end
end
%%
pha.B = B;
pha.Ay = Ay;
pha.D = D;
pha.Dind = A;
pha.bif = Dnan(A+1);
pha.numparms = {prec,nump1,nump2};
%save('Pha.mat','pha');
%%
figure(1)
imagesc(log(Ay),log(B),log(Dnan(A+1)))
%% This one
nump1 = 36;
minMu = .3;
maxMu = .7;
Mu = linspace(minMu,maxMu,nump1);

nump2 = 36;
minlogD = -6;
maxlogD = 6;
logD = linspace(minlogD,maxlogD,nump2);
D = exp(logD);


Ay = kron(10.^(-10:5),[1,2,3,5]);
Aynan = [nan,Ay]; %For disping
nump3 = length(Ay);

ax = .02;
b = 1;
c = 1;
%%
tic
A = zeros(nump2,nump1);%Index for Ay as critical value
P = zeros(nump2,nump1);%Precision utilised
TAU = zeros(nump2,nump1);%ESS
for i = 1:nump1 %loops over mu
    for j = 1:21 %loops over d
        for p = 1:5 %Loop that increases precision (and zooms) when no ESS is found
            tau = linspace(0,(b-Mu(i)),101).^p;
            ENTERED = false; %can we know this for sure?
            for k = 1:nump3
                X0 = [b-Mu(i),(b-Mu(i))/Ay(k)];
                [eq,st] = evolutionaryEq(@myModel,X0,ax,Ay(k),b,c,D(j),Mu(i),tau,1e-8);
                if ~isempty(eq)
                    TAU(i,j) = eq(end);
                    if ~ENTERED
                        [~,eqy,stability] = equilibriumsStability(ax,Ay(k),b,c,D(j),Mu(i),max(eq),true);
                        if ~isempty(stability==-1) && max(eqy(stability==-1)) > 0
                            ENTERED = true;
                        end
                    end
                elseif ENTERED
                    A(i,j) = k-1;
                    break
                end
            end
            if A(i,j) > 0 && (j == 1 || A(i,j) >= A(i,j-1))  %If A>0, then a value was found, If A(i,j) < A(i,j-1), that value is incorrect
                P(i,j) = p;
                break
            end
        end
        try
            disp(['Ay = ' num2str(Aynan(A(i,j)+1)) ' (' num2str(i) '/' num2str(nump1) ') (' num2str(j) '/' num2str(nump2) ')'])
        catch
            disp('something went wront in disp')
        end
    end
end
toc
%%
figure(2)
imagesc(log(D),Mu,log(Aynan(A+1)))
%%
pha = struct;
pha.Mu = Mu;
pha.D = D;
pha.A = A;
pha.Ay = Aynan(A+1);
pha.TAU = TAU;
save('Pha.mat','pha');