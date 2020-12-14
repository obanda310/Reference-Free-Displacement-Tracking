function [IntVal,DefVal] = CompareIntensityDeformation(directory)
if nargin ==1
    cd(directory)
end
%%
load('Matlab Data Files\FAK Adhesions.mat')
load('Matlab Data Files\3Ddata.mat','vq','image')
dSMap = double(imresize(vq.filtXY(:,:,1),size(a2f),'nearest'));
scaling = .1625; %microns per pixel
%%
mask2 = ~(bwdist(imresize(image.Area,size(a2f)))>(5/scaling)); % grab everything within 5microns of border
nMLabel = bwlabel(a2fBW);
nMLabel2 = nMLabel.*mask2;
keep = unique(nMLabel2(:));

for i = 1:max(nMLabel(:))
    if ismember(i,keep)==0
        nMLabel(nMLabel==i)=0;
    end
end

%%
dSMap2 = dSMap;
dSMap2 = dSMap2.*double(~(a2f>0));


%% 
dSMap(isnan(dSMap)) = 0;
dSMap(dSMap<.2) = 0;
dSMaxes = imregionalmax(dSMap);


%%
a3Filt = a2f>0;
SE = strel('disk',1);
a3Filt = imopen(a3Filt,SE);
a3Filt = imclose(a3Filt,SE);
a3Filt = bwareaopen(a3Filt,20);
a3Label = bwlabel(a3Filt);

%%
aProps = regionprops(nMLabel>0,img,'All');
%%
IntVal = zeros(size(aProps));
DefVal = zeros(size(aProps));
for i = 1:size(aProps,1)
    IntVal(i,1) = aProps(i).MeanIntensity;
    thisLabel = double(a3Label==a3Label(round(aProps(i).WeightedCentroid(1,2)),round(aProps(i).WeightedCentroid(1,1))));
    dSMapRegion = dSMap.*thisLabel; 
    DefVal(i,1) = max(dSMapRegion(:));
end

%%
figure
hold on
scatter(IntVal,DefVal)