A=getImages();
B=A==0;
C = bwboundaries(B);
D = zeros(size(A));
for i = 1:size(C{1,1},1)
    D(C{1,1}(i,1),C{1,1}(i,2)) = 1;
end
se = strel('disk',5');
E = uint8(255*imdilate(D,se));
imwrite(E,'BinaryOutline.tif')
figure
imshow(E)