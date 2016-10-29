clear;
close all;
clc;

% Ensures the path to necessary functions is available to the rest of the
% script
addpath(genpath('Tracking Functions'));
addpath(genpath('Kovesi Filters'));

%% Get Images and Metadata
experiment = Experiment;
images = experiment.images;
meta = experiment.metadata;
noImgs = size(images,3);
% Scaling is the same in X and Y; convert from meters to microns
pixelSize = meta.scalingX*1000000;

%% Final Pre-processing Before Finding Local Maxima
clear roiCell;
[roiImgs,roiMasks,roiCell,roiBounds,bkImg] = experiment.cropImgs;

% Create a gaussian filtered version of original to decrease false local
% maxima
ppImages7 = double(imgaussfilt(roiImgs,1.75));

% Multiply the gaussian image by the mask image to isolate regions of
% interest
ppImages8 = roiMasks.*ppImages7;

% Optional:
% ShowStack(ppImages8) %,experiment.centroids2d

%% 2D maxima approach
% Taking a break from working with 3D maxima because too many data points
% are lost in the process, and it is seeming like it will not be a good way
% to eventually identify ellipsoids and their strain/displacement. Will now
% attempt to create ellipsoids by identifying all local 2D maxima belonging
% to a single pillar, and tracing the the major axis of individual
% ellipsoids through local maxima.
clear rM cM sM rsM csM subpixMaxima tempInd tempInd2
close all

% Find local maxima in 2D (pixel resolution)
localMaxima2D = zeros(size(ppImages8));
for i = 1:noImgs
    thisMax2D = imregionalmax(ppImages8(:,:,i));    
    % In the event that a frame is empty, the local maxima are the entire
    % image (0's), this if statement removes these maxima.
    if thisMax2D == ones(size(thisMax2D,1),size(thisMax2D,2))
        thisMax2D = zeros(size(thisMax2D,1),size(thisMax2D,2));
    end
    localMaxima2D(:,:,i) = thisMax2D;
end
[rM,cM,sM] = ind2sub(size(localMaxima2D),find(localMaxima2D == 1));

% Separate maxima by z (frame) and determine indices in rM, cM, sM
for i = 1:noImgs
    if min(find(sM == i)) > 0
        twoDimMaxInd(i,1) = min(find(sM == i));
        twoDimMaxInd(i,2) = max(find(sM == i));
    end
end

% Use indices to create book of subpixel maxima for use later in linking
% maxima to pillars
for i = 1:noImgs
    if min(find(sM == i)) > 0
        [rsM,csM] = subpix2d(rM(twoDimMaxInd(i,1):twoDimMaxInd(i,2)),cM(twoDimMaxInd(i,1):twoDimMaxInd(i,2)),double(ppImages7(:,:,i)));
        subpixMaxima(1:size(rsM,2),1,i) = rsM(1,:);
        subpixMaxima(1:size(rsM,2),2,i) = csM(1,:);
        subpixMaxima(1:size(rsM,2),3,i) = i;
    end
end

%Clear out-of-bounds results from subpix2d

%-greater than x image size
[tempInd, tempInd2] = find(subpixMaxima(:,1,:) > size(ppImages8,1));
subpixMaxima(tempInd,1:3,tempInd2) = 0;
%-greater than y image size
[tempInd, tempInd2] = find(subpixMaxima(:,2,:) > size(ppImages8,2));
subpixMaxima(tempInd,1:3,tempInd2) = 0;
%-smaller than 0 in x
[tempInd, tempInd2] = find(subpixMaxima(:,1,:) < 0);
subpixMaxima(tempInd,1:3,tempInd2) = 0;
%-smaller than 0 in y
[tempInd, tempInd2] = find(subpixMaxima(:,2,:) < 0);
subpixMaxima(tempInd,1:3,tempInd2) = 0;

%% Viewing 2D pixel/subpixel maxima
% 
% close all
% %view pixel resolution maxima
% figure
% scatter3(rM,cM,sM,'.')
% 
% %view subpixel resolution maxima
% figure
% for i = 1:noImgs
%     if min(find(sM == i)) > 0
%         scatter3(subpixMaxima(:,1,i),subpixMaxima(:,2,i),subpixMaxima(:,3,i),'.')
%         hold on
%     end
% end

%% Linking objects to pillars
% The plan is to create several metrics for determining whether a local 2D
% maxima belongs to a 'pillar' group of maxima by comparing the xy distance
% between the object of interest and the nearest neighbors on frames before
% and after.
close all

% Set a maximum linking distance in microns that any object can still be
% considered part of a pillar. Smaller values will speed up code.
maxLinkDistance = 1;
maxLD = maxLinkDistance/pixelSize;

% Set a maximum number of frames to look for a linked object before giving
% up (maxJumpDistance)
maxJD = 3;

% The LinkMaxima function checks for the closest match for an object in
% later frames. Maxima with multiple matches favor pillars with a greater 
% number of constituents
[subpixMaxima,noPillars] = LinkMaxima(subpixMaxima,maxLD,maxJD);

%%
clear pillarBook
pillarSkip = 0;
skipCheck = 0;
for i = 1:noPillars
    [row,frame] = find(subpixMaxima(:,6,:)==i);
    for j = 1:size(row,1)
        if size(row,1) > 20
            pillarBook(j,1,i-pillarSkip) = subpixMaxima(row(j,1),1,frame(j,1));
            pillarBook(j,2,i-pillarSkip) = subpixMaxima(row(j,1),2,frame(j,1));
            pillarBook(j,3,i-pillarSkip) = subpixMaxima(row(j,1),3,frame(j,1));
            pillarBook(j,4,i-pillarSkip) = i-pillarSkip;
            if round(pillarBook(j,1,i-pillarSkip)) > 0 && round(pillarBook(j,2,i-pillarSkip))>0 && round(pillarBook(j,3,i-pillarSkip))>0
                pillarBook(j,5,i-pillarSkip) = roiImgs(round(pillarBook(j,1,i-pillarSkip)),round(pillarBook(j,2,i-pillarSkip)),pillarBook(j,3,i-pillarSkip));
            else
                pillarBook(j,5,i-pillarSkip) = 0;
            end
        else
            skipCheck = 1;
            break
        end
    end
    if skipCheck == 1
        pillarSkip=pillarSkip+1;
        skipCheck = 0;
    end
end
%%
createExcelForTrajectories(pillarBook);
%% 3D Scatterplot of points color coded by pillar
% figure
% for j = 1:noPillars
%     scatter3(pillarBook(:,1,j),pillarBook(:,2,j),pillarBook(:,3,j),'.')
%     hold on
% end
% hold off
%% 3D Plot of points color coded by pillar and connected
figure
for j = 1:size(pillarBook,3)
    clear tempPillar
    first = find(pillarBook(:,1,j),1,'first');
    last = find(pillarBook(:,1,j),1,'last');
    tempPillar = pillarBook(first:last,:,j);
    plot(tempPillar(:,1),tempPillar(:,2))
    hold on
end
hold off

%% Kovesi's function subpix3d
% 
% % Find local maxima in 3D (pixel resolution)
% localMaxima3D = imregionalmax(ppImages8);
% 
% % Obtain vectors with coordinates for x,y,z positions of local maxima with
% % pixel resolution
% [rM,cM,sM] = ind2sub(size(localMaxima3D),find(localMaxima3D == 1));
%
% % Find subpixel maxima based on initial guesses from imregionalmax on the
% % original images
% [rsM,csM,ssM] = subpix3d(rM,cM,sM,cropImages2);
%
% % 3D plot subpixel local maxima
% figure
% scatter3(rsM,csM,ssM,'.')

%% Dividing maxima by plane
% %Find the center of each z plane of dots based on histogram of local 3D
% %maxima
% [~,zCenters] = find(imregionalmax(imgaussfilt(histcounts(sM,noImgs),3)));
%
% % Find the spacing (# of frames) between each plane
% zSpacing = zeros(1,size(zCenters,2)-1);
% for i = 1:(size(zCenters,2)-1)
%     zSpacing(1,i) = zCenters(1,i+1)-zCenters(1,i);
% end
%
% % Choose a 'tail' size around plane centers for associating local maxima
% % above and below the plane of interest.
% zTails = round((mean(zSpacing(1,:)))/4);
%
% % Assign local maxima to planes
% zPlaneIndices = zeros(size(zCenters,2),2);
% for i = 1:size(zCenters,2)
%     zPlaneIndices(i,1) = min(find(sM>(zCenters(1,i)-zTails)));
%     zPlaneIndices(i,2) = max(find(sM<(zCenters(1,i)+zTails)));
%     zPlaneIndices(i,3) = zPlaneIndices(i,2)-zPlaneIndices(i,1);
% end
%
% zSortedMaxima = zeros(max(zPlaneIndices(:,3)),3,size(zCenters,2));
% for i = 1:size(zCenters,2)
%     zSortedMaxima(1:zPlaneIndices(i,3)+1,1,i) = rM(zPlaneIndices(i,1):zPlaneIndices(i,2),1);
%     zSortedMaxima(1:zPlaneIndices(i,3)+1,2,i) = cM(zPlaneIndices(i,1):zPlaneIndices(i,2),1);
%     zSortedMaxima(1:zPlaneIndices(i,3)+1,3,i) = sM(zPlaneIndices(i,1):zPlaneIndices(i,2),1);
% end
%
% zSortedMaximaFixed = zSortedMaxima;
% for i = 1:size(zCenters,2)
%     temp = zeros(zPlaneIndices(i,3),3);
%     temp(1:zPlaneIndices(i,3),:) = zSortedMaxima(1:zPlaneIndices(i,3),:,i);
%     fitobject{i} = fit([temp(:,1),temp(:,2)],temp(:,3),'poly11');
%     plot(fitobject{i},[temp(:,1),temp(:,2)],temp(:,3))
%     zSortedMaximaFixed(:,3,i) = fitobject{i}(zSortedMaximaFixed(:,1,i),zSortedMaximaFixed(:,1,i));
%     hold on
% end
% hold off
%
%% Predict New Local 3D maxima by Plane
%
% for i = 1:size(zCenters,2)
%     temp = zeros(zPlaneIndices(i,3),3);
%     temp(1:zPlaneIndices(i,3),:) = zSortedMaximaFixed(1:zPlaneIndices(i,3),:,i);
%     plot(fitobject{i},[temp(:,1),temp(:,2)],temp(:,3))
%     hold on
% end
% hold off
%
%%
% % Find subpixel maxima based on initial guesses from plane fits on the
% % original images
% zSortedMaximaFixedSubpix = zeros(size(zSortedMaximaFixed,1),3,size(zCenters,2));
% for i = 1:size(zCenters,2)
% [rsM2,csM2,ssM2] = subpix3d(zSortedMaximaFixed(:,1,i),zSortedMaximaFixed(:,2,i),round(zSortedMaximaFixed(:,3,i)),cropImages2);
% zSortedMaximaFixedSubpix(1:size(rsM2,2),1,i) = rsM2;
% zSortedMaximaFixedSubpix(1:size(rsM2,2),2,i) = csM2;
% zSortedMaximaFixedSubpix(1:size(rsM2,2),3,i) = ssM2;
% end
%
% for i = 1:size(zCenters,2)
%     plot(fitobject{i},[zSortedMaximaFixedSubpix(:,1,i),zSortedMaximaFixedSubpix(:,2,i)],zSortedMaximaFixedSubpix(:,3,i))
%     hold on
% end
% hold off
% figure
% for i = 1:size(zCenters,2)
% scatter3(zSortedMaximaFixedSubpix(:,1,i),zSortedMaximaFixedSubpix(:,2,i),zSortedMaximaFixedSubpix(:,3,i),'.')
% hold on
% end
% hold off