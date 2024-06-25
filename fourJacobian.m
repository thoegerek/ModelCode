function [J] = fourJacobian(SYMBOLIC)
%Jacobian for ("fourmodel")
%SYMBOLIC = bool: if function returns symbolic or function handle; false

syms n1 n2 n3 n4 a1 a2 a3 a4 b c d tau
mu = sym('mu');

f1 = b - a1*n1 -c*tau - mu;
f2 = b - a2*n2 -c*tau - mu;
f3 = b - a3*n3 -c*tau - mu;
f4 = b - a4*n4 -c*tau - mu;

L12 = -d*tau^2*n1*n2*(f1-f2);
L13 = -d*tau^2*n1*n3*(f1-f3);
L14 = -d*tau^2*n1*n4*(f1-f4);
L23 = -d*tau^2*n2*n3*(f2-f3);
L24 = -d*tau^2*n2*n4*(f2-f4);
L34 = -d*tau^2*n3*n4*(f3-f4);

d1 = f1*n1 + (f2 + mu)*n2 + (f3 + mu)*n3 + (f4 + mu)*n4 - L12 - L13 - L14;
d2 = -mu*n2 + L12 - L23 - L24;
d3 = -mu*n3 + L13 + L23 - L34;
d4 = -mu*n4 + L14 + L24 + L34;

d1d1 = diff(d1,n1);
d1d2 = diff(d1,n2);
d1d3 = diff(d1,n3);
d1d4 = diff(d1,n4);

d2d1 = diff(d2,n1);
d2d2 = diff(d2,n2);
d2d3 = diff(d2,n3);
d2d4 = diff(d2,n4);

d3d1 = diff(d3,n1);
d3d2 = diff(d3,n2);
d3d3 = diff(d3,n3);
d3d4 = diff(d3,n4);

d4d1 = diff(d4,n1);
d4d2 = diff(d4,n2);
d4d3 = diff(d4,n3);
d4d4 = diff(d4,n4);

J = [d1d1,d1d2,d1d3,d1d4;
    d2d1,d2d2,d2d3,d2d4;
    d3d1,d3d2,d3d3,d3d4;
    d4d1,d4d2,d4d3,d4d4];

if ~SYMBOLIC
    J = matlabFunction(J);
end