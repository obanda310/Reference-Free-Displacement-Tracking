function maxMasks = pillarView(images,makeMasks)
if makeMasks == 1
    roiMasks = bpass3dMB(images, [1 1 1], [7 7 11],[0 0]);
else
    roiMasks = images;
end

%% Viewing Pillars as a Frame Weighted Z-Projection
%Frames are thresholded and projected through-Z. Pixels appearing in later
%frames appear brightest.
for i = 1:size(roiMasks,3)
    temp = roiMasks(:,:,i);
    low = mean(mean(temp(temp>0)));
    
    high = max(max(roiMasks(:,:,i)));
    temp(temp>0) = high*(1);
    temp = (temp/high)*(i^4/size(roiMasks,3)^4);
    
    roiMasks2(:,:,i) = temp;
end
maxMasks = permute(max(permute(roiMasks2,[3 1 2])),[2 3 1]);
pillarView = figure;
imshow(maxMasks,[])
filePath = cd;
savefile = [filePath '\Tracking_pillarView.tif'];
export_fig(pillarView,savefile,'-native');