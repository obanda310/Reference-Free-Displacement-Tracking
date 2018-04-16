function dispSurface(directory)
if nargin ==1
    cd(directory);
end
%%
clear all
close all
tic
load('3Ddata.mat')

%%
for i = 1:size(Zeros,1)
    Zeros(i,4) = Zeros(i,1)-shear.ltLastdX(i);
    Zeros(i,5) = Zeros(i,2)-shear.ltLastdY(i);
    Zeros(i,6) = feval(Surface2,Zeros(i,4),Zeros(i,5));
end
% %%
% thenans = find(isnan(Zeros(:,1)));
% for i = 1:size(thenans,1)
%     Zeros(thenans(i),4) = shear.rawX(shear.lastFrame(thenans(i))) - shear.ltLastdX(thenans(i));
%     Zeros(thenans(i),5) = shear.rawY(shear.lastFrame(thenans(i))) - shear.ltLastdY(thenans(i));
%     Zeros(thenans(i),6) = feval(Surface2,Zeros(thenans(i),4),Zeros(thenans(i),5));
% end

%%
figure
scatter3(Zeros(:,4),Zeros(:,5),Zeros(:,6))
xlim([0 max(Zeros(:,4))+5])
ylim([0 max(Zeros(:,5))+5])
zlim([0 max(Zeros(:,6))+5])
hold on
%%
topSurface = [0,0,0,0];
surfaceFilt = image.ADil;
for i = 1:size(Zeros,1)
    %if it is under the cell
    if surfaceFilt(round((Zeros(i,2))/raw.dataKey(9,1)),round((Zeros(i,1))/raw.dataKey(9,1)))~=0
            topSurface = cat(1,topSurface,[Zeros(i,1:3) i]);        
    end
end
%shift cells up 1 to get rid of initial zero
topSurface(1,:) = [];

%%
clear dZerosZ
maxD = 1.2; % maximum positive/negative values on Z scale bar in microns
maxXY = 3; % maximum positive/negative values on XY scale bar in microns
scaleD = 32768/maxD; %scalar for creating heatmapZ
SurfCutoff = .27;

dZeros(:,1:3) = Zeros(:,1:3)-Zeros(:,4:6);
dZerosZ(:,4:6) = dZeros(:,1:3);
dZerosZ(:,1:3) = Zeros(:,4:6);
dZerosZ(topSurface(:,4),6) = 100000;
dZerosZ((abs(dZerosZ(:,6))<SurfCutoff),6) = NaN;
dZerosZ(isnan(dZerosZ(:,6)),:) = [];
dZerosZ(dZerosZ(:,6)==100000,6) = 0;

[xq,yq] = meshgrid(raw.dataKey(9,1):raw.dataKey(9,1):size(image.Black,2)*raw.dataKey(9,1), raw.dataKey(9,1):raw.dataKey(9,1):size(image.Black,1)*raw.dataKey(9,1));
vq = griddata(dZerosZ(:,1),dZerosZ(:,2),dZerosZ(:,6),xq,yq,'cubic');
xq2 = linspace(0,size(image.Black,2)*raw.dataKey(9,1),size(vq,2));
yq2 = linspace(0,size(image.Black,1)*raw.dataKey(9,1),size(vq,1));


MaximumHeatMap = imagesc(xq2,yq2,vq);
vqFinal(:,:,i) = single(vq);
imageHeat = MaximumHeatMap.CData;%.*(image.Black==0);
imageHeat(imageHeat>0) = 32768+(abs(imageHeat(imageHeat>0))*scaleD);
imageHeat(imageHeat<0) = 32768 - (abs(imageHeat(imageHeat<0))*scaleD);
imageHeat(imageHeat==0) = 32768;
imageHeat(isnan(imageHeat)) = 32768;
imageHeat = uint16(imageHeat);
imageHeatColor = single(ind2rgb(imageHeat,colorMapZ));
imageHeatColorFinal(:,:,:,i) = imageHeatColor;
close all
maxHeatMap = figure;
hold on
imshow(imageHeatColor);

%%
dZerosXY(:,1:3) = Zeros(:,4:6);
dZerosXY(:,4) = shear.ltLastdX;
dZerosXY(:,5) = shear.ltLastdY;
dZerosXY(:,6) = Zeros(:,3);

save('SurfaceData.mat','dZerosXY','dZerosZ')

