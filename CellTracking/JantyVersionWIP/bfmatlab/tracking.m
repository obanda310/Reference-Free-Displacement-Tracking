%% This needs its own section
clear;
close all;
clc;

%% Get images and metadata
[images,meta] = getCzi;
noImgs = size(images,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% OMAR ADDITIONS - 10/12/16 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Omar's preprocessing
[ppImages5,centroids2] = omarGetCentroids(images,meta);
% Clear any objects that are not 3-Dimensionally Larger than a threshold
% value
ppImages5 = bwareaopen(ppImages5,20);
ShowStack(ppImages5,centroids2)
%% Final Pre-Processing Before Finding Local Maxima
% Create ppImages 9 to fill holes in z-dimension. Double click after making
% selection to crop or an error will occur.
[ppImages6,cropImages2, rect2] = EditStack(ppImages5,images,centroids2);

%% Finding 3D Local Maxima
% Create a gaussian filtered version of original to decrease false local
% maxima
d = 3;
sig1 = 1/(1+sqrt(2))*d;
sig2 = sqrt(2) * sig1; 

% Multiply the gaussian image by the mask image to isolate regions of
% interest
ppImages8 = ppImages6.*imgaussfilt(cropImages2,sig2);
ShowStack(ppImages8,centroids2)

% Find local maxima in 3D (pixel resolution)
ppImages9 = imregionalmax(ppImages8);

%% Kovesi's Fxn subpix3d
% Obtain vectors with coordinates for x,y,z positions of local maxima with
% pixel resolution
[rM,cM,sM] = ind2sub(size(ppImages9),find(ppImages9 == 1));

% Find subpixel maxima based on initial guesses from imregionalmax on the
% original images
[rsM,csM,ssM] = subpix3d(rM,cM,sM,cropImages2);

% 3D plot subpixel local maxima
figure
scatter3(rsM,csM,ssM,'.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% END OMAR'S ADDITIONS 10/12/16 %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dividing maxima by plane
%Find the center of each z plane of dots based on histogram of local 3D
%maxima
[~,zCenters] = find(imregionalmax(imgaussfilt(histcounts(sM,noImgs),3)));

%Find the spacing (# of frames) between each plane 
zSpacing = zeros(1,size(zCenters,2)-1);
for i = 1:(size(zCenters,2)-1)
    zSpacing(1,i) = zCenters(1,i+1)-zCenters(1,i);
end

%Choose a 'tail' size around plane centers for associating local maxima
%above and below the plane of interest. 
zTails = round((mean(zSpacing(1,:)))/4);

%Assign local maxima to planes
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


%% image pre-processing
% something in here is making the features larger and blur into each other
% maybe need to get rid of noise better before doing imadjust
 imagesPro = zeros(size(images));
 for i = 1:noImgs
    imgOld = images(:,:,i);
    % Normalize image grayscale values so the image works with the
    % following built-in image processing functions.
    imgNew = imgOld/(2^meta.colorDepth-1);
    % Rolling ball background subtraction
    imgNew = rollingBall(imgNew,3);
    % Scale pixel intensities to use the entire dynamic range specified by
    % the image stack's color depth.
    imgNew = imadjust(imgNew);
    % Adjust local histograms in image to even out uneven illumination and
    % differing dot intensities within each image
    imgNew = adapthisteq(imgNew);
    % Size of structuring element radius should be JUST larger than the
    % size of the features you want to see. Dots are about 4 pixels in
    % diameter, so radius of the structuring element should be 2 or 3.
    imgNew = imtophat(imgNew,strel('disk',2));
    imagesPro(:,:,i) = imgNew;
 end