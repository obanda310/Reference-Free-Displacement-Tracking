%% Put Surface, Normal, and Empty Detections together and interpolate
%This version tries to use data from dispSurface.m
function [output] = interpFinal3D_2(directory)

if nargin ==1
    cd(directory);
end
close all

%%
load('3Ddata.mat')
%load('empties.mat') % May add this data in future versions.
load('SurfaceData.mat')
% tic
% %Create the full Lists

%This data set will try to project values from the uppermost plane to the
%surface of the hydrogel.
fullData3 = m3.ref;
fullData3(:,4:6) = m3.ref+m3.disp;
% eData = em3.ref;
% eData(:,4:6) = em3.ref+em3.disp;

%If using surfaceData
%surfaceData(:,[3 6]) = NaN;
fullData3 = cat(1,fullData3,surfaceData); %,surfaceData ,eData
%planesLoc2 = cat(2,planesLoc2,0);

%% Noise Cutoffs (Metrics from dispShear and disp3D)
% xVals = shear.ltdX(:,shear.noCellTraj);
% xVals(xVals == 0) = NaN;
% xCO = std(xVals(:),'omitnan');
% 
% yVals = shear.ltdY(:,shear.noCellTraj);
% yVals(xVals == 0) = NaN;
% yCO = std(yVals(:),'omitnan');
% 
% zCO = m3.noiseCutoff/2;
%% Rigid Body Transform (Translation followed by Rotation)
% This part of the code should make the surface of the gel be at the Z=0
% plane

%Works by creating a small plane parallel to Z=0, and then determines what
%transforms are necessary to level the plane to the to the surface approximated in
%disp3D.mat. (i.e. the hydrogel surface)

clear corner

%Create Z=0 plane
corner(1,1:2) = [0,0];
corner(2,1:2) = [1,0];
corner(3,1:2) = [0,1];
corner(4,1:2) = [1,1];

%Create gel surface parallel plane
corner(1,3) = feval(Surface2,[0,0]);
corner(2,3) = feval(Surface2,[1,0]);
corner(3,3) = feval(Surface2,[0,1]);
corner(4,3) = feval(Surface2,[1,1]);

translateZ = corner(1,3);
corner(:,3) = corner(:,3) - translateZ;

% Determine angles between planes
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

%% Translate and Rotate all of the data
disp('Zeroing coordinates to Surface')

[fD3,fD3p,fD3d] = ZeroSurfacePlane(fullData3,translateZ,tform);
for i = 1:size(plane.final,2)
    fD3(plane.final(1:nnz(plane.final(:,i)),i),4) = i;
end



fD3p = fD3(:,1:3)+fD3d;
fD3p(r.ND,:) = fD3(r.ND,1:3);
%fD3d(r.ND,:) =0;
%scatter3(fD3p(:,1),fD3p(:,2),fD3p(:,3))

fD3delete = unique(cat(1,find(isnan(fD3(:,1))),find(isnan(fD3(:,2))),find(isnan(fD3(:,3))),find(isnan(fD3d(:,1))),find(isnan(fD3d(:,2))),find(isnan(fD3d(:,3)))));
fD3(fD3delete,:) = [];
fD3p(fD3delete,:) = [];
fD3d(fD3delete,:) = [];
%fD3d(fD3(:,4)==0,3) = NaN; % Remove normal surface data (it has lower quality)
%%
figure
quiver3(fD3(:,1),fD3(:,2),fD3(:,3),fD3d(:,1),fD3d(:,2),fD3d(:,3))
%% If not using surface data, shift data up.
% clear bottomoftop
% planesLoc2(planesLoc2==0) = [];
% topPlane = find(planesLoc2 == min(planesLoc2));
% bottomoftop = min(fD3(fD3(:,4)==topPlane,3));
% if max(fD3(:,3))<0
% fD3(:,3) = fD3(:,3)-bottomoftop(1);
% end

%%
save('3Ddata.mat')
%% Newer Interp scheme 10/31/2018 (Griddata Version)
tic
dm2 = 2.12;%raw.dataKey(9,1);
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
vqZ = griddata(fD3(:,1),fD3(:,2),fD3(:,3),fD3d(:,3),xq,yq,zq);
vqZ(isnan(vqZ)) = 0;
toc

%%
  u{1}{1} = vqX * (1*10^-6);
  u{1}{2} = vqY * (1*10^-6);
  u{1}{3} = vqZ * (1*10^-6)*-1; 
  %u{1}{3}(:,:,end) = u{1}{3}(:,:,end-1);
 %%
%  u{1}{1} = vqX;
%  u{1}{2} = vqY;
%  u{1}{3} = vqZ;
%% 
save('Inputs2.mat','u')
%%
ShowStack(u{1}{1},0)
ShowStack(u{1}{2},0)
ShowStack(u{1}{3},0)
%%

dm3 = 2.12 * (1*10^-6);
  % From Example
  clear surface normals
sizeI = (size(u{1}{1}));

[surface{1}{1},surface{1}{2}] = meshgrid(dm3:dm3:sizeI(2)*dm3,dm3:dm3:sizeI(1)*dm3);
surface{1}{3} = (size(u{1}{1},3))*ones(size(surface{1}{1}))*dm3;

normals{1}{1} = zeros(size(surface{1}{1}));
normals{1}{2} = zeros(size(surface{1}{1}));
normals{1}{3} = ones(size(surface{1}{1}));


%%

model = 'linearelastic';
properties = [3900,.2];
%%
[surface, normals] = calculateSurfaceUi(surface(1), normals(1), u);
save('Inputs2.mat','u','surface','normals','model','properties','dm3')  


%%
figure
sXd = surface{2}{1}(:) -surface{1}{1}(:);
sYd = surface{2}{2}(:) -surface{1}{2}(:);
sZd = surface{2}{3}(:) -surface{1}{3}(:);
quiver3(surface{1}{1}(:),surface{1}{2}(:),surface{1}{3}(:),sXd(:),sYd(:),sZd(:),0)
%%
figure
quiver3(surface{2}{1}(:),surface{2}{2}(:),surface{2}{3}(:),normals{2}{1}(:)/100000,normals{2}{2}(:)/100000,normals{2}{3}(:)/5000000,0)
%%
try
[Fij, Sij, Eij, Uhat, ti, tiPN] = fun3DTFM(u,dm3,surface,normals,model,properties);
save('TractionOutputs.mat','Fij', 'Sij', 'Eij','Uhat','ti','tiPN')
catch
    disp('Failed Traction Coversion Function!')
end
%% Set a noise floor on Shear and Normal Tractions
tShear = imresize(tiPN{1}{1},size(image.Black));
tSFilter = tShear.*(image.ADil==0);
tSNoise = tShear.*(image.ADil~=0);
tSNoiseCO = mean(tSNoise(tSNoise(:)~=0)) + 2* std(tSNoise(tSNoise(:)~=0))
tNormal = imresize(tiPN{1}{2},size(image.Black));
tNFilter = tNormal.*(image.ADil==0);
tNNoise = tNormal.*(image.ADil~=0);
tNNoiseCO = mean(tNNoise(tNNoise(:)~=0)) + 2* std(tNNoise(tNNoise(:)~=0))
%% Noise floor v2
%%
filtI = double(imresize(image.ADil~=0,[size(tiPN{1}{1})]));
filtI(filtI==0) = NaN;
figure
imshow(tiPN{1}{1}.*filtI,[])

COShear = 2*std((tiPN{1}{1}(:).*filtI(:)),'omitnan');
CONorm = 2*std((tiPN{1}{2}(:).*filtI(:)),'omitnan');
COTotal = 2*std((ti{1}{4}(:).*filtI(:)),'omitnan');

%%
mapFilter = single(cat(3,image.ADil,image.ADil,image.ADil)==0);

mkdir('HeatMaps','Traction')
savepath = 'HeatMaps\Traction\';
figure
maxT = 175;
imshow(ti{1}{1},[])
[mapX] = Auxheatmap(size(image.Black,1),size(image.Black,2),ti{1}{1},'*blues','Xtractions',savepath,maxT,0,image.ADil);

figure
maxT = 175;
imshow(ti{1}{2},[])
[mapY] = Auxheatmap(size(image.Black,1),size(image.Black,2),ti{1}{2},'*blues','Ytractions',savepath,maxT,0,image.ADil);

figure
maxT = 300;
imshow(tiPN{1}{1},[])
[mapShear] = Auxheatmap(size(image.Black,1),size(image.Black,2),tiPN{1}{1},'*spectral','ShearTractions',savepath,maxT,COShear,image.ADil);

figure
maxT = 175;
imshow(ti{1}{3},[])
[mapZ] = Auxheatmap(size(image.Black,1),size(image.Black,2),ti{1}{3},'*blues','Ztractions',savepath,maxT,CONorm,image.ADil);

figure
maxT = 300;
imshow(tiPN{1}{2},[])
[mapNormal] = Auxheatmap(size(image.Black,1),size(image.Black,2),tiPN{1}{2},'*spectral','NormalTractions',savepath,maxT,CONorm,image.ADil);

figure
maxT = 300;
imshow(ti{1}{4},[])
[mapM] = Auxheatmap(size(image.Black,1),size(image.Black,2),ti{1}{4},'*spectral','MagnitudeTractions',savepath,maxT,COTotal,image.ADil);

%%
figure
Uhat2 = sum(Uhat{1},3);
maxT = max(Uhat2(:));
imshow(ti{1}{4},[])
Uhat2(Uhat2<1) =0 ;
[mapM] = Auxheatmap(size(image.Black,1),size(image.Black,2),Uhat2,'*spectral','StrainEnergy',savepath,maxT,0,image.ADil);

figure
Uhat3(:,:) = Uhat{1}(:,:,end);
maxT = max(Uhat3(:));
imshow(Uhat3,[])
Uhat3(Uhat3<1) =0 ;
[mapM] = Auxheatmap(size(image.Black,1),size(image.Black,2),Uhat3,'*spectral','StrainEnergyTopOnly',savepath,maxT,0,image.ADil);

%%
filtDataShear = tiPN{1}{1}.*(tiPN{1}{1}>COShear);
filtDataNorm = tiPN{1}{2}.*(tiPN{1}{2}>CONorm);
filtDataTotal = ti{1}{4}.*(ti{1}{4}>COTotal);

sumShearOld = sum(abs(ti{1}{1}(:))) + sum(abs(ti{1}{2}(:)));
sumShear = sum(abs(tiPN{1}{1}(:)));
sumShearFilt = sum(abs(filtDataShear(:)));
sumNormalOld = sum(abs(ti{1}{3}(:)));
sumNormal = sum(abs(tiPN{1}{2}(:)));
sumNormalFilt = sum(abs(filtDataNorm(:)));
sumTotal = sum(ti{1}{4}(:));
sumTotalFilt = sum(filtDataTotal(:));
NormalForce = (sumNormal*(dm3^2));
ShearForce = (sumShear*(dm3^2));
TotalForce = (sumTotal*(dm3^2));
NormalForceFilt = (sumNormalFilt*(dm3^2));
ShearForceFilt = (sumShearFilt*(dm3^2));
TotalForceFilt = (sumTotalFilt*(dm3^2));
U = sum(Uhat2(:));
Utop = sum(Uhat3(:));
output(1,1) = NormalForce;
output(1,2) = ShearForce;
output(1,3) = TotalForce;
output(1,4) = NormalForceFilt;
output(1,5) = ShearForceFilt;
output(1,6) = TotalForceFilt;
output(1,7) = U;
output(1,8) = Utop;
save('TractionStats.mat','sumShearOld','sumShear','sumNormalOld','sumNormal','U','Utop','NormalForce','ShearForce','TotalForce')

%% Write Inputs to Stackfile
for i = 1:size(u{1}{1},3)
    tempImage = uint8((u{1}{1}(:,:,i)+-1*min(min(u{1}{1}(:))))/max(max(u{1}{1}(:)))*100);
    imwrite(imresize(tempImage,size(image.Black),'nearest'),'xs.tif','WriteMode','append')
    tempImage = uint8((u{1}{2}(:,:,i)+-1*min(min(u{1}{2}(:))))/max(max(u{1}{2}(:)))*100);
    imwrite(imresize(tempImage,size(image.Black),'nearest'),'ys.tif','WriteMode','append')
    tempImage = uint8((u{1}{3}(:,:,i)+-1*min(min(u{1}{3}(:))))/max(max(u{1}{3}(:)))*100);
    imwrite(imresize(tempImage,size(image.Black),'nearest'),'zs.tif','WriteMode','append')

end

 end
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


