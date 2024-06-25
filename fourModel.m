function dX = fourModel(t,X,a1,a2,a3,a4,b,c,d,mu,tau)

n1 = max(X(1),0); n2 = X(2); n3 = X(3); n4 = X(4);

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

dX = [d1;d2;d3;d4];