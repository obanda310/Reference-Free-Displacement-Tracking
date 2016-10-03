function [centerr,centerc,intensity] = circledetectDoG(original,soluteim)
%% Detects circles using DoG filtering
% DoG filtering with additional thresholding and bqdistance algorithm to
% narrow down point selection and remove spurious detections from DoG

% Downsample rate
downsamp = 1;
im = imresize(imcomplement(original),1/downsamp);
dia = 128/downsamp;
% Calculate sigma's for DoG filter
sigma1 = 1/(1+sqrt(2)) * dia;
sigma2 = sqrt(2) * sigma1;
% Gaussfiltered images
imf1 = imgaussfilt(im,sigma1);
imf2 = imgaussfilt(im,sigma2);
% subtract filtered images
imdog = imf1-imf2;
% potenial centers are the 
allcenterpoints = imregionalmax(imdog);
 
% Limit search window
imgf = imgaussfilt(original,5);
imgfbw = im2bw(imgf,graythresh(imgf));
distim = bwdist(imgfbw);
bwd = (im2bw((distim)));
bwcircle1 = bwareafilt(bwd,2);
bwcircle2 = bwareafilt(bwcircle1,1,'smallest');
%bwdfinal = bwdist(~bwcircle2);
bwdfinal = imclose(bwcircle2,ones(30,30));
bwsmallersearchwindow = bwmorph(bwdfinal,'erode',10);
fewpoints = bwsmallersearchwindow.*allcenterpoints.*imdog;
[c,r]=find(fewpoints == max(max(fewpoints)));
centerc = round(mean(c));
centerr = round(mean(r));

%% calculate mean intensity
[imageSizeX,imageSizeY] = size(original);
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
maskradius = 60;
mask = zeros(imageSizeX,imageSizeY);
mask = (rowsInImage - centerc).^2 ...
        + (columnsInImage - centerr).^2 <= maskradius.^2;
intensity = sum(sum(soluteim.*mask));
end
