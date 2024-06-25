function out = fourLinearisedInvader(SYMBOLIC)
%The system ("fourmodel") linearised for invaders around [n1r* n2r* n3r* n4r* 0 0 0 0]
%SYMBOLIC = if function returns symbolic or function handle; false

syms n1r n2r n3r n4r n1i n2i n3i n4i a1 a2 a3 a4 b c d taur taui
mu = sym('mu');

f1 = b - a1*(n1r+n1i) -c*taui - mu;
f2 = b - a2*(n2r+n2i) -c*taui - mu;
f3 = b - a3*(n3r+n3i) -c*taui - mu;
f4 = b - a4*(n4r+n4i) -c*taui - mu;

L12 = -d*taui*n1i*(taur*n2r+taui*n2i)*(f1-f2);
L13 = -d*taui*n1i*(taur*n3r+taui*n3i)*(f1-f3);
L14 = -d*taui*n1i*(taur*n4r+taui*n4i)*(f1-f4);
L23 = -d*taui*n2i*(taur*n3r+taui*n3i)*(f2-f3);
L24 = -d*taui*n2i*(taur*n4r+taui*n4i)*(f2-f4);
L34 = -d*taui*n3i*(taur*n4r+taui*n4i)*(f3-f4);

d1 = f1*n1i + (f2 + mu)*n2i + (f3 + mu)*n3i + (f4 + mu)*n4i - L12 - L13 - L14;
d2 = -mu*n2i + L12 - L23 - L24;
d3 = -mu*n3i + L13 + L23 - L34;
d4 = -mu*n4i + L14 + L24 + L34;


d1d1 = diff(d1,n1i);
d1d2 = diff(d1,n2i);
d1d3 = diff(d1,n3i);
d1d4 = diff(d1,n4i);

d2d1 = diff(d2,n1i);
d2d2 = diff(d2,n2i);
d2d3 = diff(d2,n3i);
d2d4 = diff(d2,n4i);

d3d1 = diff(d3,n1i);
d3d2 = diff(d3,n2i);
d3d3 = diff(d3,n3i);
d3d4 = diff(d3,n4i);

d4d1 = diff(d4,n1i);
d4d2 = diff(d4,n2i);
d4d3 = diff(d4,n3i);
d4d4 = diff(d4,n4i);


o11 = subs(d1d1,{n1i,n2i,n3i,n4i},{0,0,0,0});
o12 = subs(d1d2,{n1i,n2i,n3i,n4i},{0,0,0,0});
o13 = subs(d1d3,{n1i,n2i,n3i,n4i},{0,0,0,0});
o14 = subs(d1d4,{n1i,n2i,n3i,n4i},{0,0,0,0});

o21 = subs(d2d1,{n1i,n2i,n3i,n4i},{0,0,0,0});
o22 = subs(d2d2,{n1i,n2i,n3i,n4i},{0,0,0,0});
o23 = subs(d2d3,{n1i,n2i,n3i,n4i},{0,0,0,0});
o24 = subs(d2d4,{n1i,n2i,n3i,n4i},{0,0,0,0});

o31 = subs(d3d1,{n1i,n2i,n3i,n4i},{0,0,0,0});
o32 = subs(d3d2,{n1i,n2i,n3i,n4i},{0,0,0,0});
o33 = subs(d3d3,{n1i,n2i,n3i,n4i},{0,0,0,0});
o34 = subs(d3d4,{n1i,n2i,n3i,n4i},{0,0,0,0});

o41 = subs(d4d1,{n1i,n2i,n3i,n4i},{0,0,0,0});
o42 = subs(d4d2,{n1i,n2i,n3i,n4i},{0,0,0,0});
o43 = subs(d4d3,{n1i,n2i,n3i,n4i},{0,0,0,0});
o44 = subs(d4d4,{n1i,n2i,n3i,n4i},{0,0,0,0});

J = [o11,o12,o13,o14;
    o21,o22,o23,o24;
    o31,o32,o33,o34;
    o41,o42,o43,o44];

if SYMBOLIC
    out = J;
else
    out = matlabFunction(J);
    clear n1r n2r n3r n4r n1i n2i n3i n4i a1 a2 a3 a4 b c d mu taur taui
end