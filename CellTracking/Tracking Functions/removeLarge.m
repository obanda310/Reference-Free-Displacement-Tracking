
function newMasks = removeLarge(images,masks)
subtractPattern = zeros(size(images,1),size(images,2),size(images,3));
for i = 1:size(images,3)
subtractPattern(:,:,i) = (imfill(imclose(bwmorph(edge(images(:,:,i),'Canny',[.1 .2],15),'close'),ones(50,25)),'holes')==0);
end
newMasks = subtractPattern .* masks;
end
    
