%% Reference-Free Traction Force Microscopy
% This code combines all elements of previous code into a single script.
% function RFTFM(directory)
% if nargin ==1
%     cd(directory);
% end
tic
close all
filePath = cd;
set(0,'defaultfigurecolor',[1 1 1])
addpath(genpath('Tracking Functions'));
addpath(genpath('Kovesi Filters'));
mkdir('Matlab Data Files')
%% Open Images
% First, check to see that no current image data file exists
try
    load('Matlab Data Files\CheckPoints.mat')   
catch   
    % If no image data file exists,
    auto = 1;
    disp('Creating Pre-Processed Images.')
    if auto == 0
        experiment = Experiment;
    elseif auto == 1
        [~,stackfile] = loadSpecific('*Dot Array.tif','*.tif','Did not find your file! Load a stack!');
        [~,Transfile] = loadSpecific('*Transmitted_Slice*','*.tif','Did not find your file! Load a Transmitted Image!');
        [~,Fluorfile] = loadSpecific('*FAK 576_Slice*','*.tif','Did not find your file! Load a Transmitted Image!');        
        experiment = Experiment(stackfile,Transfile,Fluorfile);
    end
    % Scaling is the same in X and Y; convert from meters to microns
    pixelSize = experiment.metadata.scalingX*1000000;
    scaleFactor = (experiment.metadata.scalingX*1000000)/0.1625;
    disp('Done.')
    disp(scaleFactor)
    disp(pixelSize)  
    
    % Process Images    
    if auto == 0
        [experiment.images,experiment.fluorImg,experiment.cellImg,roiBounds1] = experiment.cropImgs;% Initial Cropping
    end
    
    % Final Pre-processing Before Finding Local Maxima
    clear roiCell;
    [experiment.ppOptions,experiment.masks] = experiment.preprocess(auto);
    [roiImgs,roiMasks,roiCell,roiBounds,roiZeros,redoCheck] = experiment.cropImgs2(auto);
    
    while redoCheck == 1
        experiment = Experiment(experiment.images,experiment.cellImg,experiment.fluorImg,experiment.metadata);
        [experiment.ppOptions,experiment.masks] = experiment.preprocess;
        pixelSize = experiment.metadata.scalingX*1000000;
        scaleFactor = (experiment.metadata.scalingX*1000000)/0.1625;
        disp('Scale Factor:')
        disp(scaleFactor)
        [roiImgs,roiMasks,roiCell,roiBounds,roiZeros,redoCheck] = experiment.cropImgs2;
    end
    
end
%%
clear
tic
checkpoint(1,1) = 1;
save('Matlab Data Files\CheckPoints.mat','checkpoint')

% Build Image Database
autoChk = 1; %1 for auto, otherwise 0
disp('1.1 Loading Images and Raw Data')
image = ImageData(autoChk);
image = DilateBinary(image,100);
image = rawStack(image,autoChk);
image = TransDots(image);
image = maskStack(image,autoChk);
image = roiStack(image,autoChk);
image = Projection(image);

load('Matlab Data Files\DataRaw')
%Prescribed patterning settings. Can make this an input later. *Need to
%update code in several places to update to these values in order to use
%different patterning settings*
raw.pSpaceXY = 2.12;
raw.pSpaceZ = 3.5;

% Detect 3D Centroids
disp('Detecting 3D Centroids')
raw3D=RawData3D(image.MaskStack,raw);
raw3D=TranscribeR(raw3D);
disp(['done Detecting 3D Centroids at ' num2str(toc) ' seconds'])
viewDetections(raw3D,raw)
% Build Planes
disp('Building Planes')
%Set a search window size and establish neighbors
radXY =  raw.pSpaceXY*1.3; %microns
radZ = raw.pSpaceZ/9;
plane = PlanesData(raw3D);
plane = nborsPlanesF(plane,raw3D,radXY,radZ);
plane = growPlanes(plane,raw3D); %Grow from starting point until no more plane members are found
%ViewAllPlanes(plane,raw3D)%View All Planes
[plane,r] = cleanPlanes(plane,raw3D);%Filter planes with too few members (and update r)
r = RawData3D(image.MaskStack,raw,r);
r = TranscribeR(r);
ViewFilteredPlanes(plane,r)%View Filtered Planes
disp(['done Building Planes at ' num2str(toc) ' seconds'])
% Determine Likely Non-Deformed Regions
%Determine which features fall outside of a dilated mask of cell area
image = DilateBinary(image,70);
r = regionCheck(r,image.ADil,raw);
%% Build Columns
r = updatePlane(r,plane);% Flatten planes
r = updateColumn(r);% Link pillars and store info
% Build Rows
rows = RowData(r,raw,image);% Determine Row Vector
%%
r = AugmentND(r,rows);
[rows,r] = buildRows(rows,r,plane,raw,0);
disp(['done Building Rows at ' num2str(toc) ' seconds'])
[rows,r]=formatRows(rows,plane,r); % Format rows to include plane information
% Clean all major class variables
[r,rows,plane] = VarUpkeep(r,rows,plane);
%%
% View rows and verify rows were properly made
VerifyRows(r,rows)
hold on
% tidx = r.vplane==0;
% scatter3(r.X(tidx),r.Y(tidx),r.Z(tidx))
%for plotting the row Vector 'rows.V' to check alignment
pts(1,1:3) = [50+rows.V(1,1)*50 50+rows.V(2,1)*50 max(r.Z)];
pts(2,1:3) = [50-rows.V(1,1)*50 50-rows.V(2,1)*50 max(r.Z)];
plot3(pts(:,1),pts(:,2),pts(:,3))
%% Build Vertical Planes
% This part of the code should identify which rows are above/below other
% columns, and consolidate rows which have been split into multiple rows
% due to large deformations.
[VPlanes,r] = VertPlaneData(r,rows);

% Use Overlapping row data and column data to fix bad links
VPlanes = removeDuplicates(VPlanes);
% Repeat cleaning until done
clean = 1;
while clean == 1
        [VPlanes,r,rows,clean] = cleanVPlanes(VPlanes,r,rows);        
        [rows,r]=formatRows(rows,plane,r); %Update row info        
        rows = rowSizes(rows);        
        VPlanes = VertPlaneData(r,rows);
        VPlanes = removeDuplicates(VPlanes);            
end

r = updateVPlane(r,VPlanes,0);
[VPlanes,r] = ProfileColumns(VPlanes,r,plane); % Determine average pillar for reference calculations
[VPlanes,r] = AddReferenceLocs(VPlanes,r); % Add (temporary) Reference Locations to Available Candidates
[VPlanes,r] = BuildGrids(VPlanes,r,rows); % Develop grids for each VPlane
r = updateVPlane(r,VPlanes,0);
[VPlanes,r] = FillHoles(VPlanes,r,1);
disp(['done Building Vertical Planes at ' num2str(toc) ' seconds'])

% Clean all major class variables
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);
%% Verify all is working
ViewVPlanes(VPlanes,r)
%% Graphics For Publication 
ViewVPlanes3(VPlanes,r)
ViewVPlanes4(VPlanes,r)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate plane distance from surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Approximating Surface Coordinates')
[Surface2,SurfaceAll] = findSurface3(image,raw,plane);
save('Matlab Data Files\Surface.mat','Surface2')
disp(['Done Approximating Surface Coordinates at ' num2str(toc)])
try
    [plane,planesLocFiltList] = placePlanes(r,plane,Surface2);
catch
    [plane,planesLocFiltList] = placePlanes(r,plane,0);
end
%% Plot Interpolated Surface and Detections
% figure
% xlim([0 max(r.X(:))])
% ylim([0 max(r.Y(:))])
% zlim([0 size(image.RawStack,3)*raw.dataKey(10,1)])
% hold on
% for i = 1:size(plane.final,2)
%     scatter3(r.r(plane.final(1:nnz(plane.final(:,i)),i),1),r.r(plane.final(1:nnz(plane.final(:,i)),i),2),r.r(plane.final(1:nnz(plane.final(:,i)),i),3))
% end
% plot(Surface2)
% hold off
%% Improve Reference Accuracy Using Row Fits
[r,m3,image] = CalculateReferenceLocs(r,rows,plane,VPlanes,image,raw,planesLocFiltList);
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);
%% Graphics For Publication 
ViewVPlanes3(VPlanes,r)
ViewVPlanes4(VPlanes,r)
%% find average distance between markers in each row. Use those values to identify reference locations of missing objects
dXs = [];
dYs = [];
bottomRows = unique(r.row((r.planeGroup == max(r.planeGroup))));
for i = 1:size(bottomRows,1)
    thisRow = find(r.row == bottomRows(i));
    for j = 1:nnz(thisRow)
        if r.col(thisRow(j))>0
            tidx = (r.row == bottomRows(i)) & r.col ~= r.col(thisRow(j));
            [~,xdist] = dsearchn(r.rX(tidx),r.rX(thisRow(j)));
            [~,ydist] = dsearchn(r.rY(tidx),r.rY(thisRow(j)));
            if sqrt(xdist^2+ydist^2) < raw.pSpaceXY*1.5
                dXs(i,j) = xdist;
                dYs(i,j) = ydist;
            end
        end
    end
end
dXs(dXs==0) = NaN;
dYs(dYs==0) = NaN;
dXmean = mean(dXs(:),'omitnan');
dYmean = mean(dYs(:),'omitnan');

%update rowV
if rows.V(1,1) >0
    rows.V(1,1) = dXmean;
else
    rows.V(1,1) = dXmean*-1;
end


if rows.V(2,1) >0
    rows.V(2,1) = dYmean;
else
    rows.V(2,1) = dYmean*-1;
end

rows.VL = sqrt(dXmean^2+dYmean^2);

%% Build Lists of reference locations:
%
pcount = 0;
ncount = 0;
pPtsAll = [];
nPtsAll = [];
tidx = unique(r.row(r.row>0 & r.rX==0 & r.planeGroup ~=1)); %determine rows which have objects with no assigned reference position and which are not in the uppermost plane
for i = 1:nnz(tidx) %iterate through each of the identified rows
    thisRow = find(r.row==tidx(i) & r.rX>0); %find objects in the row which DO have reference positions
    for j = 1:nnz(thisRow)%for each of these objects:
        %positive direction
        pPts = rangesearch([r.rX(thisRow) r.rY(thisRow)],[r.rX(thisRow(j)) r.rY(thisRow(j))]+rows.V(1:2,1)',rows.VL/2); %determine if there is a neighbor spaced 1 row.V magnitude away
        
        %if there are none, and the object is not at the edge
        if nnz(pPts{1}) == 0 && r.rX(thisRow(j))+rows.V(1,1) > 0 && r.rX(thisRow(j))+rows.V(1,1) < r.s(2,2) && r.rY(thisRow(j))+rows.V(2,1) > 0 && r.rY(thisRow(j))+rows.V(2,1) < r.s(1,2)
            proceed = 1;
            tloc = [r.rX(thisRow(j)) r.rY(thisRow(j)) r.rZ(thisRow(j))]+rows.V(1:3,1)'; %create a reference location at the missing location
            while proceed == 1 %now search for an object to match to this reference location                
                pcount = pcount+1;                
                tidx2 = find(r.row==tidx(i) & r.rX==0); % find Markers in same row with no reference coords.
                if nnz(tidx2)> 0 %if at least 1 exists
                    [pidx pdist] = rangesearch(r.r(tidx2,1:2),tloc(1,1:2),raw.pSpaceXY/1.3); % perform a closest distance search
                    if nnz(pidx{1}) >0 % if the closest distance is within limit, assign the reference to that object
                        %pidx{1}                        
                        r.rX(tidx2(pidx{1})) = tloc(1);
                        r.rY(tidx2(pidx{1})) = tloc(2);
                        r.rZ(tidx2(pidx{1})) = tloc(3);
                    end
                    %pPtsAll = cat(1,pPtsAll,[r.rX(thisRow(j)) r.rY(thisRow(j)) r.rZ(thisRow(j))]+rows.V(1:3,1)');
                else %if none exists, break
                    proceed = 0;
                end
                %otherwise continue added row.V in the same direction to
                %check for additional missing reference locations
                tloc = tloc(1,1:3)+rows.V(1:3,1)'; % Proceed in same direction
                pPts = rangesearch([r.rX(thisRow) r.rY(thisRow)],tloc(1,1:2),rows.VL/2);
                if nnz(pPts{1}) == 0 && tloc(1)+rows.V(1) > 0 && tloc(1)+rows.V(1) < r.s(2,2) && tloc(2)+rows.V(2) > 0 && tloc(2)+rows.V(2) < r.s(1,2) % If not at edge or at other known points, proceed
                else %otherwise stop
                    proceed = 0;
                end                                
            end
        end
    end
end
r = updateColumn(r);
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);

%%
ViewVPlanes(VPlanes,r)
tidx = r.col==0 | r.row>0;
scatter3(r.rX(tidx),r.rY(tidx),r.rZ(tidx),'.')

%scatter3(pPtsAll(:,1),pPtsAll(:,2),pPtsAll(:,3))
%scatter3(nPtsAll(:,1),nPtsAll(:,2),nPtsAll(:,3))
%% Graphics For Publication 
ViewVPlanes3(VPlanes,r)
ViewVPlanes4(VPlanes,r)

%% Auto-Correct Tracking Errors
[r,rows,VPlanes,plane,~] = patchHoles(r,rows,VPlanes,plane);
[r,VPlanes,rows,plane] = clearPoorAlignment(r,VPlanes,rows,plane);
[r,rows,VPlanes,plane,trefAll] = patchHoles(r,rows,VPlanes,plane);
[VPlanes,r] = FillHoles(VPlanes,r,2);
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);

%% Graphics For Publication 
ViewVPlanes3(VPlanes,r)
ViewVPlanes4(VPlanes,r)

%% Build new scattered data to pull missing objects from
% This will use more sensitive object detection. *careful with false
% postives!
[r] = ReduceObjectDetectionThreshold(r,image,raw);
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);
%% Graphics For Publication 
ViewVPlanes3(VPlanes,r)
ViewVPlanes4(VPlanes,r)
%% Rebuild row stats
[r,m3,image] = CalculateReferenceLocs(r,rows,plane,VPlanes,image,raw,planesLocFiltList);
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);
[r,rows,VPlanes,plane,~] = patchHoles(r,rows,VPlanes,plane);
[r,VPlanes,rows,plane] = clearPoorAlignment(r,VPlanes,rows,plane);
[r,rows,VPlanes,plane,trefAll] = patchHoles(r,rows,VPlanes,plane);
[VPlanes,r] = FillHoles(VPlanes,r,2);
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);
[r,m3,image] = CalculateReferenceLocs(r,rows,plane,VPlanes,image,raw,planesLocFiltList);

% ViewRowFits(m3,r,rows,'m3')
% ViewQuiverPlot(m3,r,'m3')
%%
ViewVPlanes(VPlanes,r)
%% Graphics For Publication 
ViewVPlanes3(VPlanes,r)
ViewVPlanes4(VPlanes,r)
%% Final Error Correction
for j = 1:5
    tidx1 = find(r.planeGroup~=1 & r.colUp>0 & r.dS>0.4);
    tidx = r.colUp(tidx1((r.dS(tidx1)*1.3 > r.dS(r.colUp(tidx1))))); %incorrect matches
    tidx2 = (tidx1((r.dS(tidx1)*1.3 > r.dS(r.colUp(tidx1))))); % object below incorrect matches
    for i = 1:nnz(tidx)
        tidx3 = find(r.col==0 | r.row == 0); % Candidates to match
        tvec1 = [r.dX(tidx2(i)),r.dY(tidx2(i))];    % find appropriate vector
        tvec1(1,3) = 0;
        
        RowLineDist = zeros(size(tidx3));
        for k=1:size(tidx3,1)
            tvec2 = zeros(1,3);
            tvec2(1,1:2) = r.r(tidx3(k,1),1:2)-[r.rX(tidx(i)) r.rY(tidx(i))];
            
            RowLineDist(k,1) = norm(cross(tvec1,tvec2))/norm(tvec1); %Alignment Check 1
            RowLineDist(k,2) = dot((tvec1/norm(tvec1)),(tvec2/norm(tvec2))); %Alignment Check 2
            RowLineDist(k,3) = sqrt(sum(tvec2(1,1:2).^2)); %Distance Check 1
            RowLineDist(k,4) = sqrt(sum((r.r(tidx3(k,1),1:2)-r.r(tidx2(i),1:2)).^2));%Distance Check 2
            RowLineDist(k,5) = tidx3(k);
        end
        clear RLDCheck
        RLDCheck = RowLineDist(:,1:4);
        RLDCheck(RLDCheck(:,1)>1,1) = 0;
        RLDCheck(RLDCheck(:,2)<0.9,2) = 0;
        RLDCheck(RLDCheck(:,3)>3,3) = 0;
        RLDCheck(RLDCheck(:,4)>3,4) = 0;
        RLDCheck = double(RLDCheck>0);
        RLDCheck(:,5) = sum(RLDCheck,2);
        Match = (RowLineDist((RLDCheck(:,5)==4),5));
        
        if size(Match,1) == 1
            r.col(Match) = r.col(tidx2(i));
            r.colUp(tidx2(i)) = Match;
            r.vplane(Match) = r.vplane(tidx2(i));
            r.row(Match) = mode(VPlanes.gridRow(r.planeGroup(Match),:,r.vplane(Match)));
            r.rX(Match) = r.rX(tidx(i));
            r.rY(Match) = r.rY(tidx(i));
            r.rZ(Match) = r.rZ(tidx(i));
            r.dX(Match) = r.X(Match) - r.rX(Match);
            r.dY(Match) = r.Y(Match) - r.rY(Match);
            r.dZ(Match) = r.Z(Match) - r.rZ(Match);
            r.dS(Match) = sqrt(r.dX(Match)^2+r.dY(Match)^2);
            
            r.col(tidx(i)) = 0;
            r.row(tidx(i)) = 0;
            r.vplane(tidx(i)) = 0;
            r.rX(tidx(i)) = 0;
            r.rY(tidx(i)) = 0;
            r.rZ(tidx(i)) = 0;
            r.dX(tidx(i)) = 0;
            r.dY(tidx(i)) = 0;
            r.dZ(tidx(i)) = 0;
            r.dS(tidx(i)) = 0;
        end
        %use unpaired objects and see if any are close to the line passing
        %through the object with r.dS > the object below.
    end
end
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);
[r,rows,VPlanes,plane,~] = patchHoles(r,rows,VPlanes,plane);
[r,VPlanes,rows,plane] = clearPoorAlignment(r,VPlanes,rows,plane);
[r,rows,VPlanes,plane,trefAll] = patchHoles(r,rows,VPlanes,plane);
[VPlanes,r] = FillHoles(VPlanes,r,2);
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);
[r,m3,image] = CalculateReferenceLocs(r,rows,plane,VPlanes,image,raw,planesLocFiltList);
NoiseHists(m3,planesLocFiltList,r,'3')
ViewVPlanes(VPlanes,r)
scatter3(r.X(tidx2),r.Y(tidx2),r.Z(tidx2))
%% Graphics For Publication 
ViewVPlanes3(VPlanes,r)
ViewVPlanes4(VPlanes,r)

%% CREATE COLOR MAPS OF DISPLACEMENTS IN Z
try
  rmdir HeatMaps\3D S
catch
end
mkdir('HeatMaps\3D')
vq = DeformationMaps(r,m3,raw,plane,image);
vq = XYDeformationColorMap(vq,image,m3,raw,3,'*Spectral');
vq = ZDeformationColorMap(vq,image,m3,raw,1.2,'*PuBu');

%%
disp('Saving Data')
save('Matlab Data Files\3Ddata.mat')
disp(['Script has Completed in ' num2str(toc) ' seconds'])
%%
ViewVPlanes(VPlanes,r)
%% Unused code
%% Identify and Catalogue 2D Centroids
% % 2D dections of features should allow for identifying features which were
% % missed by the 3D object detection due to large deformations. These
% % locations will be used for locating missing features in the top plane.
% [r2] = RawData2D(image,raw);
% r2 = shape2DData(r2,raw);
%
% %% view plane 1 of 2D data
% figure
% hold on
% clear idx
% for i = 1:315
%     idx = r2.pxRaw(:,7) == i;
%     plot3(r2.pxRaw(idx,1),r2.pxRaw(idx,2),r2.pxRaw(idx,6))
% end
