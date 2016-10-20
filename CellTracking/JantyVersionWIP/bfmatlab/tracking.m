%% This needs its own section
clear;
close all;
clc;

%% Get Images and Metadata
[images,meta] = getCzi;
noImgs = size(images,3);
pixelSize = meta.scaling*1000000;

%% Pre-processing
% Process images and get centroid locations
[filtMasks,centroids] = getCentroidsStack(images,meta);
% Clear any objects that are not 3-Dimensionally larger than a threshold
% value of 20 pixels
filtMasks = bwareaopen(filtMasks,20);
ShowStack(filtMasks,centroids)

%% Final Pre-processing Before Finding Local Maxima
% Create ppImages 9 to fill hole s in z-dimension. Double click after making
% selection to crop or an error will occur.
[ppImages6,cropImages2, rect2] = EditStack(filtMasks,images,centroids);

%% Finding 3D Local Maxima
% Create a gaussian filtered version of original to decrease false local
% maxima
d = 3;
sig1 = 1/(1+sqrt(2))*d;
sig2 = sqrt(2) * sig1; 

% Multiply the gaussian image by the mask image to isolate regions of
% interest
ppImages7 = imgaussfilt(cropImages2,sig2);
ppImages8 = ppImages6.*ppImages7;
ShowStack(ppImages8,centroids)

% Find local maxima in 3D (pixel resolution)
ppImages9 = imregionalmax(ppImages8);

%% Kovesi's function subpix3d
% Obtain vectors with coordinates for x,y,z positions of local maxima with
% pixel resolution
[rM,cM,sM] = ind2sub(size(ppImages9),find(ppImages9 == 1));

% Find subpixel maxima based on initial guesses from imregionalmax on the
% original images
[rsM,csM,ssM] = subpix3d(rM,cM,sM,cropImages2);

% 3D plot subpixel local maxima
figure
scatter3(rsM,csM,ssM,'.')

%% Dividing maxima by plane
%Find the center of each z plane of dots based on histogram of local 3D
%maxima
[~,zCenters] = find(imregionalmax(imgaussfilt(histcounts(sM,noImgs),3)));

% Find the spacing (# of frames) between each plane 
zSpacing = zeros(1,size(zCenters,2)-1);
for i = 1:(size(zCenters,2)-1)
    zSpacing(1,i) = zCenters(1,i+1)-zCenters(1,i);
end

% Choose a 'tail' size around plane centers for associating local maxima
% above and below the plane of interest. 
zTails = round((mean(zSpacing(1,:)))/4);

% Assign local maxima to planes
zPlaneIndices = zeros(size(zCenters,2),2);
for i = 1:size(zCenters,2)
    zPlaneIndices(i,1) = min(find(sM>(zCenters(1,i)-zTails)));
    zPlaneIndices(i,2) = max(find(sM<(zCenters(1,i)+zTails)));
    zPlaneIndices(i,3) = zPlaneIndices(i,2)-zPlaneIndices(i,1);
end

zSortedMaxima = zeros(max(zPlaneIndices(:,3)),3,size(zCenters,2));
for i = 1:size(zCenters,2)
    zSortedMaxima(1:zPlaneIndices(i,3)+1,1,i) = rM(zPlaneIndices(i,1):zPlaneIndices(i,2),1);
    zSortedMaxima(1:zPlaneIndices(i,3)+1,2,i) = cM(zPlaneIndices(i,1):zPlaneIndices(i,2),1);
    zSortedMaxima(1:zPlaneIndices(i,3)+1,3,i) = sM(zPlaneIndices(i,1):zPlaneIndices(i,2),1);
end

zSortedMaximaFixed = zSortedMaxima;
for i = 1:size(zCenters,2)
    temp = zeros(zPlaneIndices(i,3),3);
    temp(1:zPlaneIndices(i,3),:) = zSortedMaxima(1:zPlaneIndices(i,3),:,i);
    fitobject{i} = fit([temp(:,1),temp(:,2)],temp(:,3),'poly11');
    plot(fitobject{i},[temp(:,1),temp(:,2)],temp(:,3))
    zSortedMaximaFixed(:,3,i) = fitobject{i}(zSortedMaximaFixed(:,1,i),zSortedMaximaFixed(:,1,i));
    hold on
end
hold off

%% Predict New Local 3D maxima by Plane

for i = 1:size(zCenters,2)
    temp = zeros(zPlaneIndices(i,3),3);
    temp(1:zPlaneIndices(i,3),:) = zSortedMaximaFixed(1:zPlaneIndices(i,3),:,i);
    plot(fitobject{i},[temp(:,1),temp(:,2)],temp(:,3))
    hold on
end
hold off

%%
% Find subpixel maxima based on initial guesses from plane fits on the
% original images
zSortedMaximaFixedSubpix = zeros(size(zSortedMaximaFixed,1),3,size(zCenters,2));
for i = 1:size(zCenters,2)
[rsM2,csM2,ssM2] = subpix3d(zSortedMaximaFixed(:,1,i),zSortedMaximaFixed(:,2,i),round(zSortedMaximaFixed(:,3,i)),cropImages2);
zSortedMaximaFixedSubpix(1:size(rsM2,2),1,i) = rsM2;
zSortedMaximaFixedSubpix(1:size(rsM2,2),2,i) = csM2;
zSortedMaximaFixedSubpix(1:size(rsM2,2),3,i) = ssM2;
end

for i = 1:size(zCenters,2)
    plot(fitobject{i},[zSortedMaximaFixedSubpix(:,1,i),zSortedMaximaFixedSubpix(:,2,i)],zSortedMaximaFixedSubpix(:,3,i))
    hold on
end
hold off
figure
for i = 1:size(zCenters,2)
scatter3(zSortedMaximaFixedSubpix(:,1,i),zSortedMaximaFixedSubpix(:,2,i),zSortedMaximaFixedSubpix(:,3,i),'.')
hold on
end
hold off
%% 2D maxima approach
%Taking a break from working with 3D maxima because too many data points
%are lost in the process, and it is seeming like it will not be a good way
%to eventually identify ellipsoids and their strain/displacement. Will now
%attempt to create ellipsoids by identifying all local 2D maxima belonging
%to a single pillar, and tracing the the major axis of individual
%ellipsoids through local maxima.
clear rM cM sM rsM csM ppImages10 subpixMaxima
clear tempInd tempInd2 tempInd3 tempInd4
close all
% Find local maxima in 2D (pixel resolution)
for i = 1:noImgs
    ppImages10(:,:,i) = imregionalmax(ppImages8(:,:,i));
    
    %In the event that a frame is empty, the local maxima are the entire
    %image (0's), this if statement removes these maxima.
    if ppImages10(:,:,i) == ones(size(ppImages10,1),size(ppImages10,2))
        ppImages10(:,:,i) = zeros(size(ppImages10,1),size(ppImages10,2));
    end
end
[rM,cM,sM] = ind2sub(size(ppImages10),find(ppImages10 == 1));

%Separate maxima by z (frame) and determine indices in rM, cM, sM
for i = 1:noImgs
    if min(find(sM == i)) > 0
    twoDimMaxInd(i,1) = min(find(sM == i));
    twoDimMaxInd(i,2) = max(find(sM == i));
    end
end

%use indices to create book of subpixel maxima for use later in linking
%maxima to pillars
for i = 1:noImgs 
    if min(find(sM == i)) > 0
[rsM,csM] = subpix2d(rM(twoDimMaxInd(i,1):twoDimMaxInd(i,2)),cM(twoDimMaxInd(i,1):twoDimMaxInd(i,2)),double(ppImages7(:,:,i)));
subpixMaxima(1:size(rsM,2),1,i) = rsM(1,:); 
subpixMaxima(1:size(rsM,2),2,i) = csM(1,:);
subpixMaxima(1:size(rsM,2),3,i) = i;
    end
end


%clear out-of-bounds results from subpix2d

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



%view pixel resolution maxima
figure
scatter3(rM,cM,sM,'.')

%view subpixel resolution maxima
figure
for i = 1:noImgs
    if min(find(sM == i)) > 0
    scatter3(subpixMaxima(:,1,i),subpixMaxima(:,2,i),subpixMaxima(:,3,i),'.')
    hold on
    end
end

%% Linking objects to pillars
%The plan is to create several metrics for determining whether a local 2D
%maxima belongs to a 'pillar' group of maxima by comparing the xy distance
%between the object of interest and the nearest neighbors on frames before
%and after.

%Set a maximum linking distance in microns that any object can still be
%considered part of a pillar. Smaller values will speed up code.
maxLinkDistance = 1.5;
maxLD = maxLinkDistance/pixelSize;

%Set a maximum number of frames to look for a linked object before giving
%up (maxJumpDistance)
maxJD = 2;

%Find number of objects per frame in subpixMaxima
parfor i = 1:noImgs
     [lastNZElement,~] = find(subpixMaxima(:,1,i),1,'last');
     twoDimMaxIndSub(i,1) = lastNZElement;
end


%closest neighbor in frame above
noPillars = 1
for i = 1:noImgs-1
    for j = 1:(twoDimMaxIndSub(i,1))      
        clear tempDistances
        tempDistances = sqrt((subpixMaxima(1:twoDimMaxIndSub(i,1),1,i+1)-subpixMaxima(j,1,i)).^2 +(subpixMaxima(1:twoDimMaxIndSub(i,1),2,i+1)-subpixMaxima(j,2,i)).^2);
        [nearUpNeighbor,~] = find(tempDistances==min(tempDistances));
        if tempDistances(min(nearUpNeighbor),1) < maxLD
        subpixMaxima(j,4,i) = min(nearUpNeighbor); 
        subpixMaxima(j,5,i) = tempDistances(min(nearUpNeighbor),1);
        if i>1 && ismember(j,subpixMaxima(:,4,i-1))==1
            subpixMaxima(j,6,i) = subpixMaxima(find(subpixMaxima(:,4,i-1)==j,1,'first'),6,i-1);
        else
            subpixMaxima(j,6,i) = noPillars;
            noPillars = noPillars +1
        end
        end
    end
end

% %closest neighbor in frame below
% for i = 2:noImgs
%     for j = 1:(twoDimMaxIndSub(i,1))      
%         clear tempDistances
%         tempDistances = sqrt((subpixMaxima(1:twoDimMaxIndSub(i,1),1,i-1)-subpixMaxima(j,1,i)).^2 +(subpixMaxima(1:twoDimMaxIndSub(i,1),2,i-1)-subpixMaxima(j,2,i)).^2);
%         [nearDownNeighbor,~] = find(tempDistances==min(tempDistances));
%         subpixMaxima(j,6,i) = min(nearDownNeighbor); 
%         subpixMaxima(j,7,i) = tempDistances(min(nearDownNeighbor),1); 
%     end
% end

pillarBook = zeros(noImgs,3,noPillars);
for i = 1:noPillars
    [row,frame] = find(subpixMaxima(:,6,:)==i);
    for j = 1:size(row,1)
    pillarBook(j,1,i) = subpixMaxima(row(j,1),1,frame(j,1));
    pillarBook(j,2,i) = subpixMaxima(row(j,1),2,frame(j,1));
    pillarBook(j,3,i) = subpixMaxima(row(j,1),3,frame(j,1));
    end
end

figure
for j = 500:1000
scatter3(pillarBook(:,1,j),pillarBook(:,2,j),pillarBook(:,3,j),'.')
hold on
end
hold off
