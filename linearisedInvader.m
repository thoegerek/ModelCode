function out = linearisedInvader(SYMBOLIC)
%The system ("mymodel") linearised for invaders around [xr* yr* 0 0]. Note
%the internal formulation is wrong before linearisation.
%SYMBOLIC = if function returns symbolic or function handle; false

syms xr yr xi yi ax ay b c d taur taui
mu = sym('mu');

fx = b - ax*(xr+xi) - mu - c*taui;
fy = b - ay*(yr+yi) - mu - c*taui;

D = d*taui*xi*(taur*yr + taui*yi);

dx = xi*fx + (fy + mu)*yi + D*(fx-fy);
dy = -yi*mu + D*(fy-fx); %The learning rate needs an additional term, but this should work in the limit xi,yi << 1

dxdx = diff(dx,xi);
dxdy = diff(dx,yi);
dydx = diff(dy,xi);
dydy = diff(dy,yi);

xx0 = subs(dxdx,{xi,yi},{0,0});
xy0 = subs(dxdy,{xi,yi},{0,0});
yx0 = subs(dydx,{xi,yi},{0,0});
yy0 = subs(dydy,{xi,yi},{0,0});

J = [xx0,xy0;yx0,yy0];

if SYMBOLIC
    out = J;
else
    out = matlabFunction(J);
    clear xr yr xi yi ax ay b c d mu taur taui
end