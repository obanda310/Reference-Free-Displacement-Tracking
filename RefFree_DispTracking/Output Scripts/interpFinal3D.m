%% Put Surface, Normal, and Empty Detections together and interpolate
clear all
close all

%%
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

%This data set will try to project values from the uppermost plane to the
%surface of the hydrogel.
fullData3 = m3.refSC;
fullData3(:,4:6) = m3.refSC+m3.dispFilt;

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
[fD3,fD3p,fD3d] = ZeroSurfacePlane(fullData3,translateZ,tform);
topPlane = find(planesLoc2 == min(planesLoc2));
fD3(plane.final(1:nnz(plane.final(:,topPlane)),topPlane),3) = 0;
for i = 1:size(plane.final,2)
    fD3(plane.final(1:nnz(plane.final(:,i)),i),4) = i;
    if i ~= topPlane
       fD3(plane.final(1:nnz(plane.final(:,i)),i),3) = planesLoc2(1,i)*-1; 
    end
end

fD3p = fD3(:,1:3)+fD3d;
fD3p(r.ND,:) =fD3(r.ND,1:3);
fD3d(r.ND,:) =0;
%scatter3(fD3p(:,1),fD3p(:,2),fD3p(:,3))

fD3delete = find(isnan(fD3(:,1)));
fD3(fD3delete,:) = [];
fD3p(fD3delete,:) = [];
fD3d(fD3delete,:) = [];

%% 
  clear t xqz yqz zqz fD3z
for i = 1:size(plane.final,2)
  
p1 = find(fD3(:,4) == i);
F =  scatteredInterpolant(fD3(p1,1),fD3(p1,2),fD3d(p1,3));

%Fix the grid size
gs = 2; %grid size

tx = 0:1:max(fD3p(:,1));
ty = 0:1:max(fD3p(:,2));

%Scattered X,Y to gridded x,y
[xqz{i},yqz{i}] = meshgrid(tx,ty);

%Interpolation over z-axis
z = F(xqz{i},yqz{i});

t1 = sgolayfilt(z.',4,9);
t2 = sgolayfilt(z,4,9);

t{i}  = (t1.'+t2)/2;
if i == topPlane
    zqz{i} = zeros(size(xqz{i}));
else
    zqz{i} = -1*planesLoc2(1,i)*ones(size(xqz{i}));
end
if i == 1
fD3z(:,1) = xqz{i}(:);
fD3z(:,2) = yqz{i}(:);
fD3z(:,3) = zqz{i}(:);
fD3z(:,4) = t{i}(:);
else
temp = cat(2,xqz{i}(:),yqz{i}(:),zqz{i}(:),t{i}(:));
fD3z = cat(1,fD3z,temp);
    
end

figure
surf(tx,ty,t{i})
zlim([-10,3])
colormap winter
hold on 
scatter3(fD3(p1,1),fD3(p1,2),fD3(p1,3)+fD3d(p1,3))
end

%%
figure
scatter3(fD3z(:,1),fD3z(:,2),fD3z(:,3))

%% Attempt to make all values other than surface = 0 in fD3z(:,4)
% fD3z((fD3z(:,3) ~= 0),4) = 0;
%%

%z
% fD3delete = find(abs(fD3d(:,3))<m3.noiseCutoff & abs(fD3d(:,3))~=0);
% fD3(fD3delete,:) = [];
% fD3p(fD3delete,:) = [];
% fD3d(fD3delete,:) = [];

%%
% figure
% hold on
% scatter3(fD(:,1),fD(:,2),fD(:,3))
% scatter3(fDp(:,1),fDp(:,2),fDp(:,3))

%%
% [xq,yq,zq] = meshgrid(min(fD(:,1)):raw.dataKey(9,1):max(fD(:,1)),min(fD(:,2)):raw.dataKey(9,1):max(fD(:,2)),min(fD(:,3))+2:raw.dataKey(9,1):0);
% disp('Interpolating dXs')
% vqX = griddata(fD(:,1),fD(:,2),fD(:,3),fDp(:,1)-fD(:,1),xq,yq,zq);
% vqX(isnan(vqX)) = 0;
% toc
% disp('Interpolating dYs')
% vqY = griddata(fD(:,1),fD(:,2),fD(:,3),fDp(:,2)-fD(:,2),xq,yq,zq);
% vqY(isnan(vqY)) = 0;
% toc
% disp('Interpolating dZs')
% vqZ = griddata(fDz(:,1),fDz(:,2),fDz(:,3),fDzp(:,3)-fDz(:,3),xq,yq,zq);
% vqZ(isnan(vqZ)) = 0;
% toc
%% 
tic
dm2 = raw.dataKey(9,1);
[xq,yq,zq] = meshgrid(min(fD3(:,1)):dm2:max(fD3(:,1)),min(fD3(:,2)):dm2:max(fD3(:,2)),min(fD3(:,3))+2:dm2:0);
disp('Interpolating dXs')
vqX = griddata(fD3(:,1),fD3(:,2),fD3(:,3),fD3d(:,1),xq,yq,zq);
vqX(isnan(vqX)) = 0;
toc
disp('Interpolating dYs')
vqY = griddata(fD3(:,1),fD3(:,2),fD3(:,3),fD3d(:,2),xq,yq,zq);
vqY(isnan(vqY)) = 0;
toc
disp('Interpolating dZs')
vqZ = griddata(fD3z(:,1),fD3z(:,2),fD3z(:,3),fD3z(:,4),xq,yq,zq);
vqZ(isnan(vqZ)) = 0;
toc

%% 
% ShowStack(vqX,0)
% ShowStack(vqY,0)
% ShowStack(vqZ,0)
%%
 u{1}{1} = vqX;
 u{1}{2} = vqY;
 u{1}{3} = vqZ; 
%% 
ShowStack(u{1}{1},0)
ShowStack(u{1}{2},0)
ShowStack(u{1}{3},0)
%%
dm = 1;
  %% From Example
  clear surface normals
sizeI = (size(u{1}{1})-1)*dm;

[surface{1}{1},surface{1}{2}] = meshgrid(1:dm:sizeI(2)+1,1:dm:sizeI(1)+1);
surface{1}{3} = (size(u{1}{1},3)-1)*ones(size(surface{1}{1}));

normals{1}{1} = zeros(size(surface{1}{1}));
normals{1}{2} = zeros(size(surface{1}{1}));
normals{1}{3} = ones(size(surface{1}{1}));

%%
model = 'neohookean';
properties = [12000,.2];
%%
[surface, normals] = calculateSurfaceUi(surface(1), normals(1), u);
save('Inputs2.mat','u','dm','surface','normals','model','properties')  
%%
[Fij, Sij, Eij, Uhat, ti, tiPN] = fun3DTFM(u,dm,surface,normals,model,properties);

%%
maxT = max(ti{1}{3}(:));
mkdir('HeatMaps','Traction')
savepath = 'HeatMaps\Traction\';
figure
imshow(ti{1}{1},[])
[mapX] = Auxheatmap(surface{1}{1}(1:end,1),surface{1}{2}(1,1:end),ti{1}{1},'*blues','Xtractions',savepath,maxT);

figure
imshow(ti{1}{2},[])
[mapY] = Auxheatmap(surface{1}{1}(1:end,1),surface{1}{2}(1,1:end),ti{1}{2},'*blues','Ytractions',savepath,maxT);

figure
imshow(ti{1}{3},[])
[mapZ] = Auxheatmap(surface{1}{1}(1:end,1),surface{1}{2}(1,1:end),ti{1}{3},'*blues','Ztractions',savepath,maxT);

figure
imshow(ti{1}{4},[])
[mapM] = Auxheatmap(surface{1}{1}(1:end,1),surface{1}{2}(1,1:end),ti{1}{4},'*spectral','MagnitudeTractions',savepath,maxT);


%% Functions

function [fD4,fD4p,fD4d] = ZeroSurfacePlane(fullData,translateZ,tform)
% Translate Fulldata so that the corner falls on (0,0) then rotate
% Input (:,1:6) matrix where (:,1:3) are XYZ of reference positions and
% (:,4:6) are final measured positions

% Outputs are reference position, final position, delta position each are 
% (:,1:3)
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
fD4d = fD4p - fD4;
end


