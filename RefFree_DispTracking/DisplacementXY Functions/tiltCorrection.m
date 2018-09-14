function [noiseBook,noiseStats,sumIndFinal] = tiltCorrection(roiStack,imageTrans,shear,noCellTraj)
% Tilt Correction


totalNumFrames = size(shear.rawX,1);
if nargin == 3
% sumImages = zeros(size(roiStack,1),size(roiStack,2));
% for i = 2:totalNumFrames
%     sumImages(:,:) = sumImages + roiStack(:,:,i);
% end
% sumImgScale = max(max(sumImages))/65536;
%uint16(sumImages/sumImgScale);
sumImages = uint16(squeeze(max(permute(roiStack, [3,1,2]))));
transImgScale = 65536/mean(prctile(imageTrans,95));
imageTrans = 65536-(imageTrans*transImgScale); %invert (should make opaque objects brighter)
sumImages = sumImages+imageTrans; %combine dots and cells
sumImgScale = max(max(sumImages))/65536;
sumImages = uint16(sumImages/sumImgScale);
imshow(sumImages,[]);
hold on
w = msgbox('Select a location with low displacements and double-click to continue');
                waitfor(w);
[~,sumBounds] = imcrop(sumImages);
close

sumBounds = sumBounds;
sumBounds(1,3:4) = sumBounds(1,1:2) + sumBounds(1,3:4);
sumIndX  = (shear.rawX1(:)>sumBounds(1,1) & shear.rawX1(:)<sumBounds(1,3));
sumIndY  = (shear.rawY1(:)>sumBounds(1,2) & shear.rawY1(:)<sumBounds(1,4));
sumIndXY = sumIndX .* sumIndY;
sumIndFinal = find(sumIndXY);
elseif nargin == 4
    sumIndFinal = noCellTraj;
end
disp(sumIndFinal)
shear.rawdX(shear.rawdX==0) = NaN;
shear.rawdY(shear.rawdY==0) = NaN;
shear.rawdXY(shear.rawdXY==0) = NaN;
for i = 1:totalNumFrames  
    
    noiseBook(i,1) = nanmean(shear.rawdXY(i,sumIndFinal));
    noiseBook(i,2) = nanmean(shear.rawdX(i,sumIndFinal));
    noiseStats(i,1) = nanstd(shear.rawdX(i,sumIndFinal));
    noiseBook(i,3) = nanmean(shear.rawdY(i,sumIndFinal));
    noiseStats(i,2) = nanstd(shear.rawdY(i,sumIndFinal));
    noiseBook(i,4) = sqrt(noiseBook(i,3)^2 + noiseBook(i,2)^2);
    noiseBook(i,5) = noiseBook(i,1) - noiseBook(i,4);
end

noiseBook(isnan(noiseBook)) = 0;



%if there are zeros at the top frames in noiseBook, extrapolate their
%values

x = 1:1:totalNumFrames;
weights = totalNumFrames:-1:1;
noiseFit{1} = fit(x',noiseBook(:,2),'poly1','Weights',weights');
noiseFit{2} = fit(x',noiseBook(:,3),'poly1','Weights',weights');


for i = 1:totalNumFrames
noiseStats(i,3) = feval(noiseFit{1},i)-noiseBook(i,2);
noiseStats(i,4) = feval(noiseFit{2},i)-noiseBook(i,3);
end
stdX = std(noiseStats(1:(end-5),3));
meanX = mean(noiseStats(1:(end-5),3));
stdY = std(noiseStats(1:(end-5),4));
meanY = mean(noiseStats(1:(end-5),4));
for i = 1:totalNumFrames
    if noiseStats(i,3) > (meanX+(4*stdX)) || noiseStats(i,3)<(meanX-(4*stdX))
        noiseBook(i,2) = feval(noiseFit{1},i);
    end
    if noiseStats(i,4) > (meanY+(4*stdY)) || noiseStats(i,4)<(meanY-(4*stdY))
        noiseBook(i,3) = feval(noiseFit{2},i);
    end
end



