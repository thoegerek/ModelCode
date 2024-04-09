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


D = kron(10.^(-5:5),[1,2,3,5]);
Dnan = [nan,D]; %For disping
nump3 = length(D);


%%
prec = [1, 10, 100, 1000, 10000];
A = zeros(nump1,nump2);%Index for D as critical value
P = zeros(nump1,nump2);%Precision utilised
TAU = zeros(nump1,nump2);%ESS
for i = 1:nump1 %loops over b
    for j = 1:nump2 %loops over Ay
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
                if j == 1 || TAU(i,j) <= TAU(i,j-1) %Checks for false positives arising when D large and Ay large
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