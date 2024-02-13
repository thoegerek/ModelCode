function [J,Hx,Hy] = jacobian(SYMBOLIC)
%Jacobian (and Hessians for both x and y) for ("mymodel")
%SYMBOLIC = bool: if function returns symbolic or function handle; false

syms x y ax ay b c d tau
mu = sym('mu');

fx = b - ax*x - mu - c*tau;
fy = b - ay*y - mu - c*tau;

D = d*tau^2*x*y;

dx = x*fx + (fy + mu)*y + D*(fx-fy);
dy = -y*mu + D*(fy-fx);

dxdx = diff(dx,x);
dxdy = diff(dx,y);

dydx = diff(dy,x);
dydy = diff(dy,y);

J = [dxdx,dxdy;dydx,dydy];

if ~SYMBOLIC
    J = matlabFunction(J);
end

dxdxx = diff(dxdx,x);
dxdyx = diff(dxdy,x);
dxdxy = diff(dxdx,y);
dxdyy = diff(dxdy,y);

Hx = [dxdxx,dxdxy;dxdyx,dxdyy];

dydxx = diff(dydx,x);
dydyx = diff(dydy,x);
dydxy = diff(dydx,y);
dydyy = diff(dydy,y);

Hy = [dydxx,dydxy;dydyx,dydyy];

if ~SYMBOLIC
    Hx = matlabFunction(Hx);
    Hy = matlabFunction(Hy);
    clear x y xeq yeq ax ay b c d mu tau
end