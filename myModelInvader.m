function dX = myModelInvader(t,X,ax,ay,b,c,d,mu,taur,taui)

xr = max(X(1),0); yr = X(2); 
xi = max(X(3),0); yi = X(4);

fxr = b - ax*(xr+xi) - c*taur - mu;
fyr = b - ay*(yr+yi) - c*taur - mu;
fxi = b - ax*(xr+xi) - c*taui - mu;
fyi = b - ay*(yr+yi) - c*taui - mu;

Dxr = d*taur*xr*(yr*taur + yi*taui);
Dxi = d*taui*xi*(yr*taur + yi*taui);
Dyr = d*taur*yr*(xr*taur + xi*taui);
Dyi = d*taui*yi*(xr*taur + xi*taui);

Lr = max(0,Dxr*(fyr-fxr)) - max(0,Dyr*(fxr-fyr));
Li = max(0,Dxi*(fyi-fxi)) - max(0,Dyi*(fxi-fyi));

dxr = fxr*xr + max(0,(fyr + mu)*yr) - Lr;
dyr = -mu*yr + Lr; 
dxi = fxi*xi + max(0,(fyi + mu)*yi) - Li;
dyi = -mu*yi + Li;


dX = [dxr;dyr;dxi;dyi];