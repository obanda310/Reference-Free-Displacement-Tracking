%% RUN THIS SCRIPT TO MAKE REGIONS AND OVL FILES

% OAB  April 2018
%% Adds relevant functions to path
% Works as long as folders have not been moved around
funcName = length(mfilename);
funcPath = mfilename('fullpath');
funcPath = funcPath(1:end-funcName);
cd(funcPath)
addpath(genpath([funcPath 'Dependencies']));
%% Dialogue options for all inputs
% Open images
[ListName,ListPath] = uigetfile('*.tif','Choose an image or stack to convert');
images = getImages([ListPath,ListName]);

% Choose save location and filenames. Currently set to adopt the folder
% name as the name of the files.
[TarPath] = uigetdir(ListPath,'Choose a folder to save regions to. The resulting files will inherit the folder name.');
parts = strsplit(TarPath, '\');
prefix = parts{end}; %Sets the output name prefix

% Zen Application frame parameters for scaling images and regions
pixelsize = input('What is the pixel size in microns?');
finalRes(1,1) = input('What is the final X resolution in pixels (desired frame size in ZEN)?');
finalRes(1,2) = input('What is the final Y resolution in pixels (desired frame size in ZEN)?');
tic
%% Add other options here:
%Remove small regions?
options(1,1) = 0;
%Convert single pixels to triangles
options(2,1) = 0;
%Convert single pixels to squares
options(5,1) = 1;
if options(5,1)==1
    options(2,1)=0;
end
%Run version 2 of Mask2Poly
options(3,1) = 0;
%Divide image horizontally to break open stuctures
options(4,1) = 0;

%% Iterative .Regions and .ovl outputs
for i = 1:size(images,3)
    outputName = [prefix '_' num2str(i)];
    a=images(:,:,i);
    b=a>0;
    [poly1,poly2] = Mask2Regions(b,pixelsize,outputName,TarPath,finalRes,options);
end
disp('Done!')