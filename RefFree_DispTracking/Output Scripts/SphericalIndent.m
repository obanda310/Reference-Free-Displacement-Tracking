%The goal of this code is to use displacement information from 'disp3D' to
%create interpolated displacement maps which can be directly compared to
%displacement maps generated through FE simulations.

%Coupled with an optimization script, the outputs of this code should be
%able to allow for approximation of material properties of the hydrogels
%(including the most suitable material model)using FE simulations and 
%simple indentation.

function SphericalIndent(directory)
if nargin ==1
    cd(directory);
end
%%
clear all
close all
%%
%CHANGE FONT SIZES HERE
AxisFontSize = 24;
AxisTitleFontSize = 24;
LegendFontSize = 14;

colOptions{1,1} = 'white';
colOptions{2,1} = 'black';
colOptions{1,2} = 'black';
colOptions{2,2} = 'white';
%% Making the spherical indentation data for FE fitting
load('3Ddata.mat')
%%
fullData3 = m3.ref;
fullData3(:,4:6) = m3.ref+m3.disp;
%%
 figure
 quiver3(fullData3(:,1),fullData3(:,2),fullData3(:,3),fullData3(:,4)-fullData3(:,1),fullData3(:,5)-fullData3(:,2),fullData3(:,6)-fullData3(:,3),0)
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

fD3delete = unique(cat(1,find(isnan(fD3(:,1))),find(isnan(fD3(:,2))),find(isnan(fD3(:,3))),find(isnan(fD3d(:,1))),find(isnan(fD3d(:,2))),find(isnan(fD3d(:,3)))));
fD3(fD3delete,:) = [];
fD3p(fD3delete,:) = [];
fD3d(fD3delete,:) = [];

%%
% figure
% quiver3(fD3(:,1),fD3(:,2),fD3(:,3),fD3d(:,1),fD3d(:,2),fD3d(:,3))
%% Newer Interp scheme 10/31/2018 (Griddata Version)
tic
dm2 = 2.12;%raw.dataKey(9,1);
[xq,yq,zq] = meshgrid(min(fD3(:,1)):dm2:max(fD3(:,1)),min(fD3(:,2)):dm2:max(fD3(:,2)),-20:3:-5);
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
ShowStack(vqZ)
%%
% for i =1:30
% vqtestproj = sum(vqZ,3,'omitnan')*-1;
% vqCO = prctile(vqtestproj(:),i+69);
% vqtestproj(vqtestproj<vqCO) = 0;
% 
% cmap = brewermap(25,'reds');
% 
% %inputs:
% %originalmatrix: the original matrix
% %binarisedmatrix = originalmatrix > threshold; %the thresholded matrix, a logical array
% [rows, cols] = ndgrid(1:size(vqtestproj, 1), 1:size(vqtestproj, 2));
% rc(i) = sum(rows(vqtestproj>0) .* vqtestproj(vqtestproj>0)) / sum(vqtestproj(vqtestproj>0));
% cc(i) = sum(cols(vqtestproj>0) .* vqtestproj(vqtestproj>0)) / sum(vqtestproj(vqtestproj>0));
% 
% 
% hold on
% scatter(cc(i),rc(i),50,'MarkerEdgeColor','r','linewidth',2)
% end

for i =1:5
vqtestproj = sum(vqZ,3,'omitnan')*-1;
vqCO = prctile(vqtestproj(:),i+91);
vqtestproj(vqtestproj<vqCO) = 0;

cmap = brewermap(25,'blues');

%inputs:
%originalmatrix: the original matrix
%binarisedmatrix = originalmatrix > threshold; %the thresholded matrix, a logical array
[rows, cols] = ndgrid(1:size(vqtestproj, 1), 1:size(vqtestproj, 2));
rc(i) = sum(rows(vqtestproj>0) .* vqtestproj(vqtestproj>0)) / sum(vqtestproj(vqtestproj>0));
cc(i) = sum(cols(vqtestproj>0) .* vqtestproj(vqtestproj>0)) / sum(vqtestproj(vqtestproj>0));

hold on
scatter(cc(i),rc(i),50,'MarkerEdgeColor',cmap(i,:),'linewidth',2)
end

rcf = round(mean(rc));
ccf = round(mean(cc));

%% Create map of distances to center point
distvq = zeros(size(vqZ,1),size(vqZ,2));
distvq(rcf,ccf)=1;
distvq=round(bwdist(distvq));
figure
imshow(distvq,[])
hold on
scatter(ccf,rcf,50,'MarkerEdgeColor','r','linewidth',2)

%% Create Z-Stack infographic for radial Averaging
td = max(distvq(:));

tic
for j = 1:size(vqZ,3)
    distvqmask = (abs(vqZ(:,:,j))>0).*distvq;
    for i = 1:td
        distvqA = sum(sum(distvqmask==i));
        radprofile(i,j) = sum(sum((distvq==i).*vqZ(:,:,j)))/distvqA;
    end
    toc
end

%%
    fcolor = colOptions{2,1};
    bcolor = colOptions{2,2};
set(0,'defaultfigurecolor',bcolor)
RadialProfilesPerZ=figure;
hold on
radpXX = [1:1:size(radprofile,1)]*dm2-dm2;
pgLocs = -20:3:-5;   
for i = 1:size(pgLocs,2)
    plot(radpXX,radprofile(:,i),'DisplayName',num2str(pgLocs(i)))   
end
 
    
    set(gca,'Color',bcolor,'LineWidth',2)
    %Axes, Text, Legends
    %ylim([floor(min(radprofile(:))+min(pgLocs))-1 max(pgLocs)])
    ylim([floor(min(radprofile(:)))-1 0])
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = 'Radial Distance (\mum)';% input('enter the xaxis label','s');
    yt = 'Displacement (\mum)'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    
    set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
    legend show
    legend boxoff
    leg.FontSize = LegendFontSize;
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    
    
    %Export Image
    title = ['\RadialProfilesPerZ ' fcolor ' on ' bcolor];
    savefile = [filePath title];
    export_fig(RadialProfilesPerZ,savefile,'-native');

%% Interpolate Over Set Span for Generating Heat Map.
% Try to get a 1:2 (W:L) aspect ratio
clear radprofile2 radprofile3 radprofileXs radprofileZs

radprofile2 = flip((radprofile')*-1,1);
[radprofileXs,radprofileZs] = meshgrid(radpXX,flip(pgLocs,2));
[radXX,radYY] = meshgrid(0:.2:100,[-5:-0.2:-20]);
radXX = double(radXX);
radYY = double(radYY);
radvq = griddata(radprofileXs,radprofileZs,radprofile2,radXX,radYY);

% Heatmap for XZ
figure
colorMapRad = brewermap(65536,'*spectral');
MaximumHeatMap = imagesc(radXX(:,1),radYY(:,1),radvq);
radHeat = MaximumHeatMap.CData;%.*(imageBinary==0);

radHeatNaN = (isnan(radHeat));
radHeat(isnan(radHeat)) = 0;
heatScale = (65536/(max(max(radHeat))));
radHeat = uint16(round(radHeat * heatScale));
radHeatColor = ind2rgb(radHeat,colorMapRad);

RadialHeatMap = figure;
imshow(radHeatColor)
hold on


    %Export Image
    title = ['\RadialHeatMap'];
    savefile = [filePath title];
    export_fig(RadialHeatMap,savefile,'-native');
%%
save('SphereIndent.mat','radHeatColor','radprofile2','radprofileXs','radprofileZs','radYY','radXX','radvq','cc','rc','vqtestproj','pgLocs')
disp('Spherical Indentation Script Completed Successfully')
end

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

