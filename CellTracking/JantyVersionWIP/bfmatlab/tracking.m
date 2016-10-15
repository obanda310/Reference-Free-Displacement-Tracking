%% This needs its own section
clear;
close all;
clc;

%% Get Images and Metadata
[images,meta] = getCzi;
noImgs = size(images,3);

%% Pre-processing
% Process images and get centroid locations
[filtMasks,centroids] = getCentroids(images,meta);
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
ppImages8 = ppImages6.*imgaussfilt(cropImages2,sig2);
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

for i = 1:size(zCenters,2)
    temp = zeros(zPlaneIndices(i,3),3);
    temp(1:zPlaneIndices(i,3),:) = zSortedMaxima(1:zPlaneIndices(i,3),:,i);
    fitobject{i} = fit([temp(:,1),temp(:,2)],temp(:,3),'poly11');
    plot(fitobject{i},[temp(:,1),temp(:,2)],temp(:,3))
    hold on
end

%% Predict New Local 3D maxima by Plane


