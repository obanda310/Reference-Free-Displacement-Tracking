%% Put Surface, Normal, and Empty Detections together and interpolate
clear all
close all
load('empties.mat')
load('3Ddata.mat')
load('SurfaceData.mat')
tic
%Create the full Lists


fullData = cat(1,m3.refSC,em3.refSC,dZerosXY(:,1:3));
fullData(:,4:6) = cat(1,m3.refSC+m3.dispFilt,em3.refSC+em3.dispFilt,dZerosXY(:,1:3)+dZerosXY(:,4:6));
fullData(isnan(fullData(:,1)),:) = [];

fullData2 = cat(1,m3.refSC,em3.refSC,dZerosZ(:,1:3));
fullData2(:,4:6) = cat(1,m3.refSC+m3.dispFilt,em3.refSC+em3.dispFilt,dZerosZ(:,1:3)+dZerosZ(:,4:6));
fullData2(isnan(fullData2(:,3)),:) = [];
%% Rigid Body Transform (Translation followed by Rotation)
% This part of the code should make the surface of the gel be at the Z=0
% plane

clear corner
corner(1,1:2) = [0,0];
corner(2,1:2) = [1,0];
corner(3,1:2) = [0,1];
corner(4,1:2) = [1,1];


corner(1,3) = feval(Surface2,[0,0]);
corner(2,3) = feval(Surface2,[1,0]);
corner(3,3) = feval(Surface2,[0,1]);
corner(4,3) = feval(Surface2,[1,1]);

translateZ = corner(1,3);
corner(:,3) = corner(:,3) - translateZ;
% Generate Angles
Ya = -atan(corner(2,3));
Xa = atan(corner(3,3));
Za = 0;
EulerA = zeros(4,4);
EulerA(1,1) = cos(Ya);
EulerA(1,2) = sin(Za)*sin(Ya);
EulerA(1,3) = sin(Ya) *cos(Za);
EulerA(2,1) = sin(Xa) *sin(Ya);
EulerA(2,2) = cos(Xa)*cos(Za)-cos(Ya)*sin(Xa)*sin(Za);
EulerA(2,3) = -cos(Xa)*sin(Za)-cos(Ya)*cos(Za)*sin(Xa);
EulerA(3,1) = -cos(Xa)*sin(Ya);
EulerA(3,2) = cos(Za)*sin(Xa)+cos(Xa)*cos(Ya)*sin(Za);
EulerA(3,3) = cos(Xa)*cos(Ya)*cos(Za)-sin(Xa)*sin(Za);
EulerA(4,4) = 1;

% Test Euler Proper Angles
corner2 = pointCloud(corner);
tform = affine3d(EulerA);
corner3 = pctransform(corner2,tform);

%%
disp('Zeroing coordinates to Surface')
[fD,fDp] = ZeroSurfacePlane(fullData,translateZ,tform);
[fDz,fDzp] = ZeroSurfacePlane(fullData2,translateZ,tform);

%%
% figure
% hold on
% scatter3(fD(:,1),fD(:,2),fD(:,3))
% scatter3(fDp(:,1),fDp(:,2),fDp(:,3))

%%
[xq,yq,zq] = meshgrid(min(fD(:,1)):raw.dataKey(9,1):max(fD(:,1)),min(fD(:,2)):raw.dataKey(9,1):max(fD(:,2)),min(fD(:,3))+2:raw.dataKey(9,1):0);
disp('Interpolating dXs')
vqX = griddata(fD(:,1),fD(:,2),fD(:,3),fDp(:,1)-fD(:,1),xq,yq,zq);
vqX(isnan(vqX)) = 0;
toc
disp('Interpolating dYs')
vqY = griddata(fD(:,1),fD(:,2),fD(:,3),fDp(:,2)-fD(:,2),xq,yq,zq);
vqY(isnan(vqY)) = 0;
toc
disp('Interpolating dZs')
vqZ = griddata(fDz(:,1),fDz(:,2),fDz(:,3),fDzp(:,3)-fDz(:,3),xq,yq,zq);
vqZ(isnan(vqZ)) = 0;
toc
%% 
ShowStack(vqX,0)
ShowStack(vqY,0)
ShowStack(vqZ,0)
%%

[xq2,yq2] = meshgrid(1:dm:size(u{1}{1},2)+1,1:dm:size(u{1}{1},1)+1);
zq2 = ones(size(xq2))*35;

% u{1}{1} = permute(vqX,[2,1,3]);
% u{1}{2} = permute(vqY,[2,1,3]);
% u{1}{3} = permute(vqZ,[2,1,3]);
%%
u{1}{1} = permute(u{1}{1},[2,1,3]);
u{1}{2} = permute(u{1}{2},[2,1,3]);
u{1}{3} = permute(u{1}{3},[2,1,3]);

dm = 1;

  %% From Example
  clear surface normals
sizeI = (size(u{1}{1})-1)*dm;

[surface{1}{1},surface{1}{2}] = meshgrid(1:dm:sizeI(2)+1,1:dm:sizeI(1)+1);
surface{1}{3} = 18*ones(size(surface{1}{1}));

normals{1}{1} = zeros(size(surface{1}{1}));
normals{1}{2} = zeros(size(surface{1}{1}));
normals{1}{3} = ones(size(surface{1}{1}));


model = 'linearElastic';
properties = [12000,.35];

[surface, normals] = calculateSurfaceUi(surface(1), normals(1), u);
save('Inputs2.mat','u','dm','surface','normals','model','properties')  
%%
[Fij, Sij, Eij, Uhat, ti, tiPN] = fun3DTFM(u,dm,surface,normals,model,properties);

%%
figure
imshow(tiPN{1}{1},[])

figure
imshow(tiPN{1}{2},[])

%% Functions

function [fD4,fD4p] = ZeroSurfacePlane(fullData,translateZ,tform)
% Translate Fulldata so that the corner falls on (0,0) then rotate
fD2 = fullData(:,1:3);
fD2(:,3) = fD2(:,3)-translateZ;
fD2 = pointCloud(fD2(:,1:3));

fD2p(:,1:3) = fullData(:,4:6);
fD2p(:,3) = fD2p(:,3)-translateZ;
fD2p = pointCloud(fD2p(:,1:3));
%rotate
fD3 = pctransform(fD2,tform);
fD3p = pctransform(fD2p,tform);

fD4 = fD3.Location;
fD4p = fD3p.Location;
end


