clear images imagesFilled imagesHoles imagesBigHoles imagesSmallHoles imagesD imagesWS imagesFilled2
close all
se = strel('disk',3);
se2 = strel('disk',3);

images=  imerode(imdilate(imerode(imdilate(bpass(image3,3,50)>0,se),se),se),se2);

imagesFilled = imfill(images,'holes');
imagesHoles = imagesFilled & ~images;
for i = 1:25
    imagesBigHoles(:,:,i) = bwareaopen(imagesHoles,i*40);
    imagesSmallHoles(:,:,i) = imagesHoles & ~imagesBigHoles(:,:,i);
    imagesD(:,:,i) = -1*bwdist(~(images | imagesSmallHoles(:,:,i)),'quasi-euclidean');
    imagesD(~(images | imagesSmallHoles(:,:,i)))= Inf;
    imagesWS(:,:,i) = double(watershed(imagesD(:,:,i))>0);
    imagesFilled2(:,:,i) = double((images | imagesSmallHoles(:,:,i)));
end
ShowStack(imagesFilled2)
ShowStack(imagesWS.*imagesFilled2)