function out = linearised(SYMBOLIC)
%The system ("mymodel") linearised for base population
%SYMBOLIC = if function returns symbolic or function handle; false

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

if SYMBOLIC
    out = J;
else
    out = matlabFunction(J);
    clear x y ax ay b c d mu tau
end