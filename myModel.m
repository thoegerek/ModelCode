function dX = myModel(t,X,ax,ay,b,c,d,mu,tau)

x = max(X(1),0); y = X(2);

fx = b - ax*x -c*tau - mu;
fy = b - ay*y -c*tau - mu;

D = d*x*y*tau^2;

dx = fx*x + max((fy + mu)*y,0) + D*(fx-fy);
dy = -mu*y + D*(fy-fx); 

dX = [dx;dy];