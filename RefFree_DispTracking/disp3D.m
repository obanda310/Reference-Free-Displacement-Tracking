function disp3D(directory)
if nargin ==1
    cd(directory);
end
clear all
close all
set(0,'defaultfigurecolor',[1 1 1])
tic
autoChk = 1; %1 for auto, 0 otherwise
%% Load Shear Data
load('Shear Mat Files\DataRaw')
load('Shear Mat Files\DataShear')

%% Load Images
image = ImageData(autoChk);
image = rawStack(image,autoChk);
image.RawStack=permute(image.RawStack, [2,1,3]);
image = TransDots(image);
%imshow(image.TDots,[]);

%% Kilfoil Stack Filter
disp('Processing Images')
res=bpass3dMB(image.RawStack, [1 1 1], [7 7 11],[0 0]);
%ShowStack(res)
disp(['done Processing Images at '  num2str(toc) ' seconds'])

%% Kilfoil Object Detection 3D

disp('Detecting 3D Centroids')
raw3D=RawData3D(res,raw);
raw3D=TranscribeR(raw3D);
disp(['done Detecting 3D Centroids at ' num2str(toc) ' seconds'])
viewDetections(raw3D,raw)

% Build Planes Dot by Dot
%%
disp('Building Planes')
%Set a search window size and establish neighbors
radXY =  2.5; %microns
radZ = .3;
plane = PlanesData(raw3D);
plane = nborsPlanesF(plane,raw3D,radXY,radZ);
%Grow from starting point until no more plane members are found
plane = growPlanes(plane,raw3D);

disp(['done Building Planes at ' num2str(toc) ' seconds'])

%%
% View all detected planes
figure
hold on
for i = 1:size(plane.raw,2)
    scatter3(raw3D.X(plane.raw(1:nnz(plane.raw(:,i)),i)),raw3D.Y(plane.raw(1:nnz(plane.raw(:,i)),i)),raw3D.Z(plane.raw(1:nnz(plane.raw(:,i)),i)))
end
hold off
%%
% Filter planes with too few members (and update r)
[plane,r] = cleanPlanes(plane,raw3D);
r = RawData3D(res,raw,r);
r = TranscribeR(r);
toc
%%
%View Filtered Planes
pf = figure;
hold on
for i = 1:size(plane.final,2)
    scatter3(r.X(plane.final(1:nnz(plane.final(:,i)),i)),r.Y(plane.final(1:nnz(plane.final(:,i)),i)),r.Z(plane.final(1:nnz(plane.final(:,i)),i)))
end
fcolor = 'white';
bcolor = 'black';
AxisFontSize = 24;
LegendFontSize = 14;
xt = 'X';% input('enter the xaxis label','s');
    yt = 'Y'; %input('enter the yaxis label','s');
    zt = 'Z';
    label{1} = xlabel(xt);
    label{2} = ylabel(yt);
    label{3} = zlabel(zt);
    set(gca,'YMinorTick','on','color',bcolor)
    ytickformat('%.1f')
        le{1} = 'plane 1'; %input('enter the legend','s');
    le{2} = 'plane 2';
    ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,1,[0 0])
    %errorbar(meanDisplacements(1,1:3),meanDisplacements(2,1:3),'.','color',[0 0 0],'MarkerSize',1)
    axis([0 max(r.X) 0 max(r.Y) 0 ceil(max(r.Z))])
    legend off
hold off

 savefile = 'XZ Indent.tif';
    export_fig(pf,savefile,'-native');
%% Dots in Cell Region
%Determine which features fall outside of a dilated mask of cell area
image = DilateBinary(image,50);
r = regionCheck(r,image.ADil,raw);

%% Check to see if a previous row slope exists
clear files check
files = dir('*.mat'); %Check Directory for default filenames
if size(files,1)>=1
    for k = 1:size(files,1)
        current=files(k).name;
        check(k)=strcmp(current(end-7:end),'rowV.mat');
    end
    loc=find(check);
else
    loc= zeros(1,1);
end
disp('Generating Row Slope')
image = FindNDSquare(image);
if size(loc,2)>0
    load('rowV.mat')
    disp('Found a previous row slope!')
else
    [rowV,rowV2, rROI, rROIplane] = rowVselector(r,raw,image);
    %rowV = rowV2;
end
rowV(1,3) = 0;
disp(['done Generating Row Slope at ' num2str(toc) ' seconds'])

%% Assign detections to Rows
disp('Building Rows')
[r,rows] = buildRows2(r,rowV,plane.final);
disp(['done Building Rows at ' num2str(toc) ' seconds'])
% figure
% hold on
% for i = 1:size(rows,1)
%     scatter3(r.X(rows(i,1:nnz(rows(i,:)))),r.Y(rows(i,1:nnz(rows(i,:)))),r.Z(rows(i,1:nnz(rows(i,:)))))
% end

%% Format rows to include plane information
[rows,rowPlanes,rowPlanesIdx,rowsNDCU,r] = formatRows(rows,plane,r);
% figure
% hold on
% for j = 1%1:size(rowPlanes,3)
%     for i = 1:size(rowPlanes,1)
%         n = nnz(rowPlanes(i,:,j));
%         scatter3(r.X(rowPlanes(i,1:n,j)),r.Y(rowPlanes(i,1:n,j)),r.Z(rowPlanes(i,1:n,j)))
%     end
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate plane distance from surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Approximating Surface Coordinates')
[Surface2,SurfaceAll,Zeros] = findSurface2(shear,image,raw);
disp(['Done Approximating Surface Coordinates at ' num2str(toc)])
try
    [planesLoc2,planesLocFiltList] = placePlanes(r,plane,Surface2);
catch
    [planesLoc2,planesLocFiltList] = placePlanes(r,plane,0);
end
%% Plot Interpolated Surface and Detections
% figure
% xlim([0 max(shear.rawX(:))])
% ylim([0 max(shear.rawY(:))])
% zlim([0 size(image.RawStack,3)*raw.dataKey(10,1)])
% hold on
% for i = 1:size(plane.final,2)
%     scatter3(r.r(plane.final(1:nnz(plane.final(:,i)),i),1),r.r(plane.final(1:nnz(plane.final(:,i)),i),2),r.r(plane.final(1:nnz(plane.final(:,i)),i),3))
% end
% plot(SurfaceAll)
% hold off
%% Match 3D detections to 2D-based pillars
[r] = match2D(r,raw,shear,rowsNDCU);
% figure
% hold on
% for i=1:max(r.col(:))
%     current = find(r.col(:)==i);
%     if size(current,1)>0
%         plot3(r.X(current),r.Y(current),r.Z(current))
%     end
%
% end
% for i = 1:size(shear.rawX,2)
% scatter3(shear.rawX(:,i),shear.rawY(:,i),shear.rawZ(:,i),'.')
% end
%
%%
%
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
%ViewMethodRef(m1,r,'m1')
m1.disp = r.r(:,1:3)-m1.refSC(:,1:3);
% Displacement Statistics for Planes at least 4.5 microns from surface
% Find all dots farther than 4.5 microns from surface
% Calculate Average Displacement per Row per Plane Method 1
m1 = dispStats(m1,plane,rowPlanesIdx,r,rowsNDCU,planesLocFiltList);
% Filter out noise in Displacements Method 1
m1 = dispNoise(m1,r,planesLocFiltList,image,shear,raw.dataKey(9,1));
disp(['done Fitting Rows with Independent Slopes - Method 1 at ' num2str(toc) ' seconds'])
% Scatter3/Plot3 of Dots/Fits Method 2
%ViewRowFits(m1,r,rowPlanes,'m1')
% Quiver Plot of Displacements Method 2
%ViewQuiverPlot(m1,r,'m1')
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
%ViewMethodRef(m2,r,'m2')
m2.disp = r.r(:,1:3)-m2.refSC(:,1:3);
% Calculate Average Displacement per Row per Plane Method 2
m2 = dispStats(m2,plane,rowPlanesIdx,r,rowsNDCU,planesLocFiltList);
% Filter out noise in Displacements Method 1
m2 = dispNoise(m2,r,planesLocFiltList,image,shear,raw.dataKey(9,1));
disp(['done Fitting Rows with Row Slope - Method 2 at ' num2str(toc) ' seconds'])
% Scatter3/Plot3 of Dots/Fits Method 2
%ViewRowFits(m2,r,rowPlanes,'m2')
% Quiver Plot
%ViewQuiverPlot(m2,r,'m2')
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
m3.disp = r.r(:,1:3)-m3.refSC(:,1:3);
% Calculate Average Displacement per Row per Plane Method 3
m3 = dispStats(m3,plane,rowPlanesIdx,r,rowsNDCU,planesLocFiltList);
% Filter out noise in Displacements Method 1

m3 = dispNoise(m3,r,planesLocFiltList,image,shear,raw.dataKey(9,1));
image.imgNBds = imageNoiseBounds(r.r,m3,image,shear,raw.dataKey(9,1),m3.noiseCutoff);
disp(['done Picking Best Case Row Fit - Method 3 at ' num2str(toc) ' seconds'])
%%
%ViewMethodRef(m3,r,'m3')
%Scatter3/Plot3 of Dots/Fits Method 3
%ViewRowFits(m3,r,rowPlanes,rowPlanesIdx,'m3')
%Quiver Plot of Displacements Method 3
%ViewQuiverPlot(m3,r,'m3')
%%
NoiseHists(m3,planesLocFiltList,r,'3')


% %%
% figure
% axis([0 r.s(1,2) 0 r.s(2,2)])
% hold on
% for j = 1:size(rowPlanes,3)
%     if nnz(plane.raw(:,j))>5000
%         m3fit{j} = fit([r(plane.raw(1:nnz(plane.raw(:,j)),j),1),(size(image.Area,1)*raw.dataKey(9,1))-r(plane.raw(1:nnz(plane.raw(:,j)),j),2)],r(plane.raw(1:nnz(plane.raw(:,j)),j),3),'lowess','Span',.02);
%         plot(m3fit{j})
%     end
%     plot(fitSurface)
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE COLOR MAPS OF DISPLACEMENTS IN Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir('HeatMaps\3D')
mkdir('HeatMaps\3D\ColorBar')
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
savefile = [filePath '\HeatMaps\3D\ColorBar\ColorBarZ.tif'];
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
savefile = [filePath '\HeatMaps\3D\ColorBar\ColorBarShear.tif'];
export_fig(colorBarSave,savefile,'-native');

SE = strel('disk',round(10/.1625));
imageBinaryCombined = image.imgNBds; %imdilate(image.ADil==0,SE);%+(imageBinary2==0))==0;
%
%Determine which 'planes' in plane.final should be the same plane
clear planesGroups
for i = 1:size(planesLoc2,2)
    clear differences
    differences = planesLoc2 - planesLoc2(1,i);
    planesGroups(i,1:size(find(abs(differences)<2),2)) = find(abs(differences)<2)';
end
planesGroups = unique(planesGroups,'rows');
[HeatMapN,vqN] = heatmapZ(r.r,m3.disp,plane.final,planesGroups,imageBinaryCombined,image.Borders,raw.dataKey(9,1),0,colorMapZ,colorMapXY,0);
%[HeatMap,vq1] = heatmapZ(r.r,m1.disp,plane.final,planesGroups,imageBinaryCombined,image.Borders,raw.dataKey(9,1),m1.noiseCutoff,colorMapZ,colorMapXY,1);
%[HeatMap2,vq2] = heatmapZ(r.r,m2.disp,plane.final,planesGroups,imageBinaryCombined,image.Borders,raw.dataKey(9,1),m2.noiseCutoff,colorMapZ,colorMapXY,2);
[HeatMap3,vq3] = heatmapZ(r.r,m3.dispFilt,plane.final,planesGroups,imageBinaryCombined,image.Borders,raw.dataKey(9,1),m3.noiseCutoff,colorMapZ,colorMapXY,3);
[HeatMap3,vq4] = heatmapZ(r.r,m3.disp,plane.final,planesGroups,imageBinaryCombined,image.Borders,raw.dataKey(9,1),m3.noiseCutoff,colorMapZ,colorMapXY,32);


% Print Data to txt file
% Calculate Useful parameters
%[vq1pos,vq1neg,m1.PosTotal, m1.PosMax, m1.NegTotal, m1.NegMax]  = vqStats(vq1,planesGroups);
%[vq2pos,vq2neg,m2.PosTotal, m2.PosMax, m2.NegTotal, m2.NegMax]  = vqStats(vq2,planesGroups);
[vq3pos,vq3neg,m3.PosTotal, m3.PosMax, m3.NegTotal, m3.NegMax]  = vqStats(vq3,planesGroups);

% Write to file
cd HeatMaps
planesLocTxt = fopen('Average Z location of planes.txt','wt');
p1Format = 'Plane no. %1.0f is at %.2f microns from the surface \n';
%p2Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 1 fit \n';
%p3Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 2 fit \n';
p4Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 3 fit \n';
for i = 1:size(planesGroups,1)
    fprintf(planesLocTxt,p1Format,i,mean(planesLoc2(1,planesGroups(i,(planesGroups(i,:)>0)))));
    %fprintf(planesLocTxt,p2Format,i,mean(m1.MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))));
    %fprintf(planesLocTxt,p3Format,i,mean(m2.MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))));
    fprintf(planesLocTxt,p4Format,i,mean(m3.MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))));
end
fclose(planesLocTxt);
cd(filePath)

%Comma Delimiter Version
cd HeatMaps
%printStats(m1,planesGroups,planesLocTxt,planesLoc2,'1');
%printStats(m2,planesGroups,planesLocTxt,planesLoc2,'2');
printStats(m3,planesGroups,planesLocTxt,planesLoc2,'3');
cd(filePath)


%%
disp('Saving Data')
save 3Ddata
%%

% Save data for profile views
%%
filePath = cd;
folderName = 'Profile Data';
mkdir(filePath,folderName)
save('Profile Data\vqZ.mat','vqN','vq3','image','HeatMapN','HeatMap3')

%%
save('rowV.mat','rowV')

%% Attempt to interpolate normal displacements in 3D
% rDispFilt = rDisp;
% rDispFilt((rDispFilt<.4 & rDispFilt>-.4))=0;
%
% meshRes = 2;
% %[xx,yy,zz] = meshgrid(1:meshRes:size(image.RawStack,1),1:meshRes:size(image.RawStack,2),1:1:size(image.RawStack,3));
% [xx,yy,zz] = meshgrid(1:meshRes:r.s(1),1:meshRes:r.s(2),min(r.Z(:)):0.2:max(r.Z(:)));
% xq = double([xx(:) yy(:) zz(:)]);
% vq = griddatan(rRef(:,1:3),rDispFilt(:,3),xq);
% vq = reshape(vq,size(xx));
% vq(isnan(vq)) = min(min(min(vq)));
%
% vq2 = vq;
% vq2 = vq2+abs(min(min(min(vq2))));
% vq2Scale = double(65000/max(max(max(vq2))));
% vq2 = uint16(vq2Scale*vq2);
% ShowStack(vq2,1,1)

disp(['Script has Completed in ' num2str(toc) ' seconds'])