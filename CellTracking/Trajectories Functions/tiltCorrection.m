function [noiseBook,sumIndFinal] = tiltCorrection(roiStack,book1,book2,dataKey)
% Tilt Correction
close all
totalNumFrames = size(book1,2);
sumImages = zeros(size(roiStack,1),size(roiStack,2));
for i = 2:totalNumFrames
    sumImages(:,:) = sumImages + roiStack(:,:,i);
end
sumImgScale = max(max(sumImages))/65536;
sumImages = uint16(sumImages/sumImgScale);
imshow(sumImages,[]);
hold on
w = msgbox('Select a location with low displacements and double-click to continue');
                waitfor(w);
[~,sumBounds] = imcrop(sumImages);
close

sumBounds = sumBounds*dataKey(7,1);
sumBounds(1,3:4) = sumBounds(1,1:2) + sumBounds(1,3:4);
sumIndX  = (book2(:,1)>sumBounds(1,1) & book2(:,1)<sumBounds(1,3));
sumIndY  = (book2(:,2)>sumBounds(1,2) & book2(:,2)<sumBounds(1,4));
sumIndXY = sumIndX .* sumIndY;
sumIndFinal = find(sumIndXY);

%book1(book1==0) = NaN;
for i = 1:totalNumFrames  
    noiseBook(i,1) = nanmean(book1(5,i,sumIndFinal));
    noiseBook(i,2) = nanmean(book1(3,i,sumIndFinal));
    noiseBook(i,3) = nanmean(book1(4,i,sumIndFinal));
    noiseBook(i,4) = sqrt(noiseBook(i,3)^2 + noiseBook(i,2)^2);
    noiseBook(i,5) = noiseBook(i,1) - noiseBook(i,4);
end
%book1(isnan(book1)) = 0;