function dX = myModelInvader(t,X,ax,ay,b,c,d,mu,taur,taui)

xr = max(X(1),0); yr = X(2); 
xi = max(X(3),0); yi = X(4);

fxr = b - ax*(xr+xi) - c*taur - mu;
fyr = b - ay*(yr+yi) - c*taur - mu;
fxi = b - ax*(xr+xi) - c*taui - mu;
fyi = b - ay*(yr+yi) - c*taui - mu;

Dr = d*taur*xr*(yr*taur + yi*taui);
Di = d*taui*xi*(yr*taur + yi*taui);

dxr = fxr*xr + (fyr + mu)*yr + Dr*(fxr-fyr);
dyr = -mu*yr + Dr*(fyr-fxr); 
dxi = fxi*xi + (fyi + mu)*yi + Di*(fxi-fyi);
dyi = -mu*yi + Di*(fyi-fxi);


dX = [dxr;dyr;dxi;dyi];