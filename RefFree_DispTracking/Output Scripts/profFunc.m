function [normXY,normZ,normAxis,cell_boundary] = profFunc(directory)
if nargin == 1
cd(directory)
end
clear all 
close all
%%
set(0,'defaultfigurecolor',[1 1 1])
load('Profile Data\vqXY.mat');
load('Profile Data\HeatMapXY.mat');
load('Profile Data\vqZ.mat');
load('3Ddata.mat','planesLoc2','planesGroups')
zTarget = 7;
for j = 1:size(planesGroups,1)
        planesLoc3(j) = mean(planesLoc2(1,planesGroups(j,1:nnz(planesGroups(j,:)))));
end
zPlane = find(abs(planesLoc3-zTarget) == min(abs(planesLoc3-zTarget)),1,'first');
vqXY = imresize(vqXY,[size(vq3,1) size(vq3,2)]);
imageHeatXYColor = imresize(imageHeatXYColor,[size(vq3,1) size(vq3,2)]);
%% open bw area image
files = dir('*.tif'); %Check Directory for default filenames
filePath = strcat(cd,'\');
clear check
    for k = 1:length(files)
    current=files(k).name;
    if length(current)>=15 
    check(k)=strcmp(current(end-14:end),'Binary Mask.tif');
    end
    end
    loc=find(check);
    if size(loc,1)==1
    imageArea= imread(files(loc(1)).name);
    imageArea = imresize(imageArea,(size(image.Trans,1)/size(imageArea,1)));
    else
    [nameAreaFile,filePath] = uigetfile('*.tif','Select a Thresholded Image of the Cell Area');
    imageArea = imread([filePath,nameAreaFile]);
    imageArea = imresize(imageArea,(size(image.Trans,1)/size(imageArea,1)));
    end

imageAreaPos = double(imageArea==0);
imageAreaPosFilt = bwareaopen(imageAreaPos,5000);
imACentroid = regionprops(imageAreaPosFilt,'Centroid');
[xMax,yMax] = find(vqXY == max(max(vqXY)));
xCent = round(imACentroid.Centroid(1,1));
yCent = round(imACentroid.Centroid(1,2));
%%
[normXY,normZ,normAxis,dispProfXY,dispProfZ,heatMap,cell_boundary] = profileDisp(HeatMapN(:,:,:,zPlane),imageHeatXYColor,vqXY,vqN(:,:,zPlane),[xCent yCent],[yMax xMax],imageArea,HeatMap3(:,:,:,1));
