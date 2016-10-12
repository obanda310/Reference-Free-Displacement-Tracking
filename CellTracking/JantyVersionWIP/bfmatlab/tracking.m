%%
clear all
close all
%% import image stack
% choose .czi image stack
[name,path] = uigetfile('*.czi');
filename = [path,name];
% read into matlab
czi = bfopen(filename);
% get image dimensions from first image in .czi file
[imgRows, imgCols] = size(czi{1}{1,1});
% get number of images (z-slices) from .czi file
global noImgs 
noImgs= length(czi{1});
% column 1 in cell 1 contains all images, row index corresponds to position
% in z stack
images = zeros(imgRows,imgCols,noImgs);
parfor i = 1:noImgs
    images(:,:,i) = czi{1}{i,1};
end
%% read metadata
metadata = czi{1,2};
binning          = metadata.get('Global Information|Image|Channel|Binning #1');
zStart           = metadata.get('Global Experiment|AcquisitionBlock|MultiTrackSetup|ZStackSetup|First|Distance|Value #1'); %um
zEnd             = metadata.get('Global Experiment|AcquisitionBlock|MultiTrackSetup|ZStackSetup|Last|Distance|Value #1'); %um
exposureTime     = metadata.get('Global HardwareSetting|ParameterCollection|ExposureTime #1');
colorDepth       = metadata.get('Global Information|Image|ComponentBitCount #1');
scaling          = metadata.get('Global Scaling|Distance|Value #1');
global sizeX;sizeX = str2num(metadata.get('Global Information|Image|SizeX #1'));
global sizeY;sizeY = str2num(metadata.get('Global Information|Image|SizeY #1'));

%% omar's preprocessing
close all
global centroids2
centroids2= cell(noImgs,1);
d=3;
sig1 = 1/(1+sqrt(2))*d;
sig2 = sqrt(2) * sig1;
noiseRng = mean(prctile(images(:,:,noImgs),95));
[mVal,mIdx] = max(prctile(mean(images(:,:,:)),70));
refHistImage = images(:,:,mIdx);
LoG = fspecial('log',3,.25);
Lap = fspecial('laplacian');
% imadjust(ppImages, [noiseRng;max(max(max(images)))],[0,1];
parfor i = 1:noImgs
highIn = max(max(max(images(:,:,i))));
ppImages(:,:,i) = (imadjust(images(:,:,i)/65535,[noiseRng/65535 highIn/65535],[])*65535);
ppImages2(:,:,i) = (imgaussfilt(imadjust(images(:,:,i)/65535,[noiseRng/65535 highIn/65535],[])*65535,sig2)-imgaussfilt(imadjust(images(:,:,i)/65535,[noiseRng/65535 highIn/65535],[])*65535,sig1))*10;
ppImages3(:,:,i) = imfilter(ppImages2(:,:,i),Lap);
end
ppImages4 = (ppImages3>(prctile(ppImages3(ppImages3>0),25)));


parfor i = 1:noImgs
ppImages5(:,:,i) = bwareaopen(ppImages4(:,:,i),5);
    c = regionprops(ppImages5(:,:,i),'Centroid');
    centroids2{i} = cat(1,c.Centroid);
end

%Clear any objects that are not 3-Dimensionally Larger than a threshold
%value
ppImages5 = bwareaopen(ppImages5,20);
ShowStack(ppImages5)
%% Final Pre-Processing Before Finding Local Maxima
%create ppImages 9 to fill holes in z-dimension
ppImages6 = ppImages5;

%% Finding 3D Local Maxima
%create a gaussian filtered version of original to decrease false local
%maxima
parfor i = 1:noImgs
    ppImages7(:,:,i) = imgaussfilt(images(:,:,i),sig2);
end

%multiply the gaussian image by the mask image to isolate regions of
%interest
ppImages8 = ppImages6.*ppImages7;

ShowStack(ppImages8)

%find local maxima in 3D (pixel resolution)
ppImages9 = imregionalmax(ppImages8);

%% Kovesi's Fxn subpix3d
%obtain vectors with coordinates for x,y,z positions of local maxima with
%pixel resolution
[rM,cM,sM] = ind2sub(size(ppImages9),find(ppImages9 == 1));

%find subpixel maxima based on initial guesses from imregionalmax on the
%original images
[rsM,csM,ssM] = subpix3d(rM,cM,sM,images);

%3D plot subpixel local maxima
figure
scatter3(rsM,csM,ssM,'.')

%% image pre-processing
% something in here is making the features larger and blur into each other
% maybe need to get rid of noise better before doing imadjust
 imagesPro = zeros(size(images));
 parfor i = 1:noImgs
    imgOld = images(:,:,i);
    % Normalize image grayscale values so the image works with the
    % following built-in image processing functions.
    imgNew = imgOld/(2^str2double(colorDepth)-1);
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
%% paper steps
% 1. define point spread function
%   a. empirically
%   b. analytically
% 2. focal adhesion detection
%   a. background subtraction
%   b. contrast limited adaptive histogram equalization
%   c. laplacian of gaussian filter
%   d. manual thresholding
% 3. dot detection
%   a. xy
%       i. threshold
%       ii. dot postion calculted via weighted centroid of grayscale values
%           in feature
%   b. z
%       i. spot detection function of Imaris
%       ii. fit a plane through all points, take difference between plane
%           and measured z coordinates = tilt correction
% 4. mesh generation (some other time)