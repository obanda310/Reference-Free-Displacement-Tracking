function [subtractPattern] = removeLarge(images) %newMasks,
subtractPattern = zeros(size(images,1),size(images,2),size(images,3));
SE = strel('disk',3,0);
for i = 1:size(images,3)
    subtractPattern(:,:,i) = imerode(imfill(imdilate(bwmorph(edge(images(:,:,i),'Canny',[.1 .2],15),'close'),SE),'holes'),SE)==0;
end
%newMasks = subtractPattern .* masks;
end
    
