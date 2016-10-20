[images,meta] = getCzi;
noImgs = size(images,3);
imagesNew = zeros(size(images));
parfor i = 1:noImgs
    thisImg = images(:,:,i);
    masks(:,:,i) = createMask(thisImg);
%     imagesNew(:,:,i) = thisImg .* masks(:,:,i);
end
%%
centroidDots = imregionalmax(imagesNew);
[rM,cM,sM] = ind2sub(size(centroidDots),find(centroidDots == 1));
%%
% Find subpixel maxima based on initial guesses from imregionalmax on the
% original images
[rsM,csM,ssM] = subpix3d(rM,cM,sM,imagesNew);
%%
% 3D plot subpixel local maxima
figure
scatter3(rsM,csM,ssM,'.')