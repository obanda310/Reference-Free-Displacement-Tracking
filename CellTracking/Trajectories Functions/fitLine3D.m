function [A,B] = fitLine3D(x,y,z)

n=size(x,1);

%3D line fit
P=[mean(x),mean(y),mean(z)]';
[~,~,V]=svd([x-P(1),y-P(2),z-P(3)]);
N=1/V(end,1)*V(:,1);


A=P+dot([x(1),y(1),z(1)]'-P,N)*N/norm(N)^2;
B=P+dot([x(n),y(n),z(n)]'-P,N)*N/norm(N)^2;