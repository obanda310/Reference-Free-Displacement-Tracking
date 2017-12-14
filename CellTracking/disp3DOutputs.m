function disp3DOutputs(directory,redo)
tic
close all
if nargin == 1 || nargin == 2
cd(directory);
end
clear all
load('3Ddata.mat')
%% Redo Displacement calculations
if nargin == 2 && redo ==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st METHOD FOR MEASURING DISPLACEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Fitting Rows with Independent Slopes - Method 1')
% Extend the fit Lines passed the edge of the dot area Method 1
m1 = DispData3D;
m1 = method1fit(m1,r,rows,rowV);
% Calculate Displacements from Fit Lines Method 1
m1 = calcDisp(m1,r,rows,rowsNDCU); %based on fits of non-deformed markers
m1 = calcDispSC(m1,r,shear); %includes shift correction in reference approximation
%ViewMethodRef(m1,r)
m1.disp = r.r(:,1:3)-m1.refSC(:,1:3);
% Displacement Statistics for Planes at least 4.5 microns from surface
% Find all dots farther than 4.5 microns from surface
% Calculate Average Displacement per Row per Plane Method 1
m1 = dispStats(m1,plane,rowPlanesIdx,r,rowsNDCU,planesLocFiltList);
% Filter out noise in Displacements Method 1
m1 = dispNoise(m1,r,planesLocFiltList);
disp('done Fitting Rows with Independent Slopes - Method 1')
% Scatter3/Plot3 of Dots/Fits Method 2
%ViewRowFits(m1,r,rowPlanes)
% Quiver Plot of Displacements Method 2
%ViewQuiverPlot(m1,r)
NoiseHists(m1,planesLocFiltList,r,'1')



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd METHOD FOR MEASURING DISPLACEMENTS V2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Fitting Rows with Row Slope - Method 2')
m2 = DispData3D;
m2 = method2fit(m2,m1,r,raw,image.ADil,rows,rowPlanesIdx);
% Calculate Displacements from Fit Lines Method 2
m2 = calcDisp(m2,r,rows,rowsNDCU); %based on fits of non-deformed markers
m2 = calcDispSC(m2,r,shear); %includes shift correction in reference approximation
%ViewMethodRef(m2,r)
m2.disp = r.r(:,1:3)-m2.refSC(:,1:3);
% Calculate Average Displacement per Row per Plane Method 2
m2 = dispStats(m2,plane,rowPlanesIdx,r,rowsNDCU,planesLocFiltList);
% Filter out noise in Displacements Method 1
m2 = dispNoise(m2,r,planesLocFiltList);
disp('done Fitting Rows with Row Slope - Method 2')
% Scatter3/Plot3 of Dots/Fits Method 2
%ViewRowFits(m2,r,rowPlanes)
% Quiver Plot
%ViewQuiverPlot(m2,r)
NoiseHists(m2,planesLocFiltList,r,'2')



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3RD METHOD FOR FITTING ROWS (METHOD 1 AND 2 COMBINED)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Picking Best Case Row Fit - Method 3')
m3 = DispData3D;
m3 = method3fit(m3,m1,m2,rows);
% Calculate Displacements from Fit Lines Method 3
m3 = calcDisp(m3,r,rows,rowsNDCU); %based on fits of non-deformed markers
m3 = calcDispSC(m3,r,shear); %includes shift correction in reference approximation
%ViewMethodRef(m3,r)
m3.disp = r.r(:,1:3)-m3.refSC(:,1:3);
% Calculate Average Displacement per Row per Plane Method 3
m3 = dispStats(m3,plane,rowPlanesIdx,r,rowsNDCU,planesLocFiltList);
% Filter out noise in Displacements Method 1
m3 = dispNoise(m3,r,planesLocFiltList);
disp('done Picking Best Case Row Fit - Method 3')
% Scatter3/Plot3 of Dots/Fits Method 3
%ViewRowFits(m3,r,rowPlanes)
% Quiver Plot of Displacements Method 3
%ViewQuiverPlot(m3,r)
NoiseHists(m3,planesLocFiltList,r,'3')

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE COLOR MAPS OF DISPLACEMENTS IN Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorMapZ = single(brewermap(65536,'*PuBu'));
colorBar1 = single(zeros(500,25));
range = uint16(round(linspace(65536,1,500)'));
for i = 1:25
    colorBar1(1:500,i) = range;
end
colorBar2 = ind2rgb(colorBar1,colorMapZ);
for i = 1:10
    colorBar2((i*50)-3:(i*50),13:25,:) = 0;
end
colorBar2(1:3,13:25,:) = 0;

%Save Color Bar Image
close all
colorBarSave = figure;
hold on
imshow(colorBar2);
filePath=cd;
savefile = [filePath '\HeatMaps\ColorBarZv2.tif'];
export_fig(colorBarSave,savefile,'-native');

colorMapXY = single(brewermap(65536,'*Spectral'));
colorBar1 = single(zeros(500,25));
range = uint16(round(linspace(65536,1,500)'));
for i = 1:25
    colorBar1(1:500,i) = range;
end
colorBar2 = single(ind2rgb(colorBar1,colorMapXY));
for i = 1:10
    colorBar2((i*50)-3:(i*50),13:25,:) = 0;
end
colorBar2(1:3,13:25,:) = 0;

%Save Color Bar Image
close all
colorBarSave = figure;
hold on
imshow(colorBar2);
filePath=cd;
savefile = [filePath '\HeatMaps\ColorBarZv2.tif'];
export_fig(colorBarSave,savefile,'-native');

SE = strel('disk',round(10/.1625));
imageBinaryCombined = imdilate(image.ADil==0,SE);%+(imageBinary2==0))==0;
%%
%Determine which 'planes' in plane.final should be the same plane
clear planesGroups
for i = 1:size(planesLoc2,2)
    clear differences
    differences = planesLoc2 - planesLoc2(1,i);
    planesGroups(i,1:size(find(abs(differences)<2),2)) = find(abs(differences)<2)';
end
planesGroups = unique(planesGroups,'rows');
[HeatMapN,vqN] = heatmapZ(r.r,m3.disp,plane.final,planesGroups,imageBinaryCombined,raw.dataKey(9,1),0,colorMapZ,colorMapXY,0);
[HeatMap3,vq3] = heatmapZ(r.r,m3.dispFilt,plane.final,planesGroups,imageBinaryCombined,raw.dataKey(9,1),m3.noiseCutoff,colorMapZ,colorMapXY,3);

% Print Data to txt file
% Calculate Useful parameters
[vq3pos,vq3neg,m3.PosTotal, m3.PosMax, m3.NegTotal, m3.NegMax]  = vqStats(vq3,planesGroups);

% Write to file
cd HeatMaps
planesLocTxt = fopen('Average Z location of planes.txt','wt');
p1Format = 'Plane no. %1.0f is at %.2f microns from the surface \n';
p4Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 3 fit \n';
for i = 1:size(planesGroups,1)
    fprintf(planesLocTxt,p1Format,i,mean(planesLoc2(1,planesGroups(i,(planesGroups(i,:)>0)))));
    fprintf(planesLocTxt,p4Format,i,mean(m3.MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))));
end
fclose(planesLocTxt);
cd(filePath)

%Comma Delimiter Version
cd HeatMaps
printStats(m3,planesGroups,planesLocTxt,planesLoc2,'3');
cd(filePath)

%%
save 3Ddata
disp(['Script has Completed in ' num2str(toc) ' seconds'])
% Save data for profile views
%%
filePath = cd;
folderName = 'Profile Data';
mkdir(filePath,folderName)
save('Profile Data\vqZ.mat','vqN','vq1','vq2','vq3','image','HeatMap','HeatMap2','HeatMap3')