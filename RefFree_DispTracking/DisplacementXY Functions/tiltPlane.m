function [tiltPlaneFitMicrons,tiltPlaneFitPixels] = tiltPlane(noiseBook,dataKey,imageArea)% test data

%Uses tilt information to calculate the photopatterning plane relative to
%the imaging plane (an alternative to fitting ellipsoid centers)

clear n P N x y z
n=size(noiseBook,1);
x=(noiseBook(:,2)+(size(imageArea,1)/2));
y=(noiseBook(:,3)+(size(imageArea,2)/2));
z=(linspace(1,n,n)')*dataKey(10,1);

%3D line fit
P=[mean(x),mean(y),mean(z)]';
[~,~,V]=svd([x-P(1),y-P(2),z-P(3)]);
N=1/V(end,1)*V(:,1);


A=P+dot([x(1),y(1),z(1)]'-P,N)*N/norm(N)^2;
B=P+dot([x(n),y(n),z(n)]'-P,N)*N/norm(N)^2;

%three points on perpendicular plane for a fit
C = 1000*null((B-A).');
C1(:,1) = C(:,1)+B;
C1(:,2) = B;
C1(:,3) = C(:,2)+B;

tiltPlaneFitMicrons{1} = fit([C1(2,:)',C1(1,:)'],C1(3,:)','poly11');

C1(1:2,:) = C1(1:2,:)/dataKey(9,1);
C1(3,:) = C1(3,:)/.4;

tiltPlaneFitPixels{1} = fit([C1(2,:)',C1(1,:)'],C1(3,:)','poly11');
