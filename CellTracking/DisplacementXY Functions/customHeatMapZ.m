function [imageHeat ,vq] = customHeatMapZ(shear,imageBlack,dataKey,outputs,filePath)
disp('7.4 Creating Surface Normal Deformation Heat Map')
%Create Color Bar
colorMap2 = brewermap(65536,'*RdGy');
colorBar1 = zeros(500,25);
range = uint16(round(linspace(65536,1,500)'));
for i = 1:25
    colorBar1(1:500,i) = range;
end
colorBar2 = ind2rgb(colorBar1,colorMap2);
for i = 1:10
    colorBar2((i*50)-3:(i*50),13:25,:) = 0;
end
colorBar2(1:3,13:25,:) = 0;

%Save Color Bar Image
close all
colorBarSave = figure;
hold on
imshow(colorBar2);
maxD = 1.5; % maximum positive/negative values on scale bar in microns
scaleD = 32768/maxD; %scalar for creating heatmap
savefile = [filePath '\HeatMaps\ColorBarZ.tif'];
if ismember(6,outputs) == 1
    export_fig(colorBarSave,savefile,'-native');
else
    export_fig(colorBarSave,savefile);
end
close

%Create Heat Map of Max Displacement
[xq,yq] = meshgrid(0:2*dataKey(9,1):size(imageBlack,2)*dataKey(9,1), 0:2*dataKey(9,1):size(imageBlack,1)*dataKey(9,1));
vq =  griddata(shear.rawX1(:),shear.rawY1(:),shear.dTop2(:),xq,yq,'cubic');
xq2 = linspace(0,size(imageBlack,2)*dataKey(9,1),size(vq,2));
yq2 = linspace(0,size(imageBlack,1)*dataKey(9,1),size(vq,1));
MaximumHeatMap = imagesc(xq2,yq2,vq);
imageHeat = MaximumHeatMap.CData;
imageHeat(imageHeat>0) = 32768+(abs(imageHeat(imageHeat>0))*scaleD);
imageHeat(imageHeat<0) = 32768 - (abs(imageHeat(imageHeat<0))*scaleD);
imageHeat(isnan(imageHeat)) = 32768;
imageHeat = uint16(imageHeat);

imageHeatColor = ind2rgb(imageHeat,colorMap2);

%Save Max Heat Map Image
close all
maxHeatMap = figure;
hold on
imshow(imageHeatColor);
savefile = [filePath '\HeatMaps\Single\MaximumHeatMapZ.tif'];
if ismember(6,outputs) == 1
    export_fig(maxHeatMap,savefile,'-native');
else
    export_fig(maxHeatMap,savefile);
end
close

%Save Color Bar Values (depends on vq!)
cd HeatMaps
colorBarTxt = fopen('Heat Map Color Bar Values Z.txt','wt');
p1Format = 'Tick number %1.0f is %.2f \n';
for i = 1:11
    colorBarValues(i,1) =  round(maxD - ((maxD/5)*(11-i)),2); % round(max(vq(:)) - ((max(vq(:))/10)*(11-i)),2)
    fprintf(colorBarTxt,p1Format,i,colorBarValues(i,1));
end
fclose(colorBarTxt);
cd(filePath)
end