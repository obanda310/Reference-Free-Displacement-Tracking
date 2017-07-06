clear all
close all
image = getImages();

%% Create Thresholded Image
images2 = zeros(size(image,1),size(image,2));
for i = 1:50
    image2 = bpass(image,3,i);
    images2 = images2 +image2;
end
image6 = images2/50;
image7 = image6<prctile(mean(image6),99);

%% Fill small holes
se = strel('disk',9);
image8 = imerode(image7,se);
imagesBigHoles = bwareaopen(image8==0,1500);
imagesSmallHoles = image8==0 & ~imagesBigHoles;
image9 = image7|imagesSmallHoles;

%% Identify Regions
se2 = strel('disk',2);
image9e = imerode(image9,se2);
image10 = -1*bwdist(~image9e);
image10g = imgaussfilt(image10,7);
image10i = image10g;
image10i(~(image9e | imagesSmallHoles))= Inf;
imagesWS = double(watershed(image10i)>0);
image11 = image9.*imagesWS;
image12 = bwlabel(image11);
imshow(image12,[])
%% Remove Regions on Edges of Image
EdgeSize = 50;
image13 = image12;
for i = 1:max(max(image13))
    temp = image13 == i;
    if sum(sum(temp(1:EdgeSize,:))) > 0 || sum(sum(temp(:,1:EdgeSize))) > 0 || sum(sum(temp((size(temp,1)-EdgeSize):size(temp,1),:))) > 0 || sum(sum(temp(:,(size(temp,2)-EdgeSize):size(temp,2)))) > 0 
        image13(image13==i) = 0;
    end
end
image13 = bwlabel(image13>0);
imshow(image13,[])

%% Obtain statistics on remaining regions
clear regionsStats
regionStats = regionprops(image13>0,'Area','Centroid','Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength','Perimeter');
%% View regions
clear regions
for i = 1:max(max(image13))
    regions(:,:,i) = image13 == i;
end
ShowStack(regions)
%% Histograms!
bins = 20;

figure
subplot(2,2,1)
histogram(cat(1,regionStats.EquivDiameter),bins)
title('Diameter, pixels')
subplot(2,2,2)
histogram(cat(1,regionStats.Area),bins)
title('Area, pixels')
subplot(2,2,3)
histogram(cat(1,regionStats.Perimeter),bins)
title('Perimeter, pixels')
subplot(2,2,4)
histogram(cat(1,regionStats.MajorAxisLength),bins)
title('Major Axis, pixels')

