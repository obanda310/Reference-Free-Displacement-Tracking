function [imageHeat2max,vq2max,vq2total,imageHeatNaN,imageHeatColorM]=customHeatMap(shear,deformation,image,dataKey,outputs,filePath)

%Input Variables to Specify
maxD = 3; %maximum value on the heat map color bar (microns)
res = 1; % number of original pixels between interpolated points (value of 1 should return a heat map with original image dimensions)
res2 = res;

% Create Edge Clearing Mask
Edges = (image.ADil == 0);

% Creating Heat Maps!
disp('7.1 Creating Maximum Shear Heat Maps')
close all
clear xq yq vq xq2 yq2
folderName = ('HeatMaps');
%cleardir = rmdir(folderName,'s');
mkdir(filePath,folderName)
folderName = ('Shear');
mkdir(strcat(filePath,'\HeatMaps\'),folderName)

%Create Color Bar
colorMap2 = brewermap(65536,'*spectral');
%colorMap2(1,1:3) = 0;
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
savefile = [filePath '\HeatMaps\Shear\MaximumColorBar.tif'];
if ismember(6,outputs) == 1
    export_fig(colorBarSave,savefile,'-native');
else
    export_fig(colorBarSave,savefile);
end
close

%Create Heat Map of Max Displacement
[xq,yq] = meshgrid(1:res:size(Edges,2), 1:res:size(Edges,1));
xq = xq*dataKey(9,1);
yq = yq*dataKey(9,1);
vq = dataKey(9,1)*griddata(shear.rawX1(:),shear.rawY1(:),deformation(:)/dataKey(9,1),xq,yq,'linear');
%vq = vq.*Edges;
vq2max = griddata(shear.rawX1(:),shear.rawY1(:),shear.coFilt(:),xq,yq,'linear');
vq2max = vq2max.*Edges;
disp(num2str(max(max(vq))))
xq2 = linspace(0,size(Edges,2)*dataKey(9,1),size(vq,2));
yq2 = linspace(0,size(Edges,1)*dataKey(9,1),size(vq,1));

MaximumHeatMap = imagesc(xq2,yq2,vq);
imageHeat = MaximumHeatMap.CData;
imageHeatNaN = (isnan(imageHeat));
imageHeat(isnan(imageHeat)) = 0;
heatScale = (65536/maxD); %(max(max(imageHeat)))
imageHeat = uint16(round(imageHeat * heatScale));
imageHeatColor = ind2rgb(imageHeat,colorMap2);


MaximumHeatMap2 = imagesc(xq2,yq2,vq2max);
imageHeat2 = MaximumHeatMap2.CData;
imageHeat(isnan(imageHeat)) = 0;
imageHeat2 = uint16(round(imageHeat2 * heatScale));
imageHeatColor2 = ind2rgb(imageHeat2,colorMap2);
imageHeat2max = imageHeat2;
imageHeatColorM = imageHeatColor2;

%Cutout the Area of cell or pattern
MaximumHeatMap3 = imagesc(xq2,yq2,vq);
vq3 = double(vq).*double(image.Area>0);
imageHeat3 = MaximumHeatMap3.CData;
imageHeat3(isnan(imageHeat3)) = 0;
imageHeat3 = uint16(round(imageHeat3 * heatScale));
imageHeat3 = imageHeat3.*uint16(image.Area>0);
imageHeatColor3 = ind2rgb(imageHeat3,colorMap2);
imageHeatColor3M = imageHeatColor3;

save([filePath '\Shear Mat Files\ShearCutout.mat'],'vq3')


%Save Max Heat Map Image
close all
maxHeatMap = figure;
hold on
imshow(imageHeatColor);

maxHeatMap2 = figure;
hold on
imshow(imageHeatColor2);

maxHeatMap3 = figure;
hold on
imshow(imageHeatColor3);

if ismember(6,outputs) == 1
    savefile = [filePath '\HeatMaps\Shear\MaximumHeatMap.tif'];
    export_fig(maxHeatMap,savefile,'-native');
    savefile = [filePath '\HeatMaps\Shear\MaximumHeatMapNoiseCutoff.tif'];
    export_fig(maxHeatMap2,savefile,'-native');
        savefile = [filePath '\HeatMaps\Shear\MaximumHeatMapBinaryCutout.tif'];
    export_fig(maxHeatMap3,savefile,'-native');
else
    savefile = [filePath '\HeatMaps\Shear\MaximumHeatMap.tif'];
    export_fig(maxHeatMap,savefile);
    savefile = [filePath '\HeatMaps\Shear\MaximumHeatMapNoiseCutoff.tif'];
    export_fig(maxHeatMap2,savefile);
    savefile = [filePath '\HeatMaps\Shear\MaximumHeatMapBinaryCutout.tif'];
    export_fig(maxHeatMap3,savefile);
end




%Save Color Bar Values (depends on vq!)
cd HeatMaps\Shear
colorBarTxt = fopen('Heat Map Color Bar Values.txt','wt');
p1Format = 'Tick number %1.0f is %.2f \n';
for i = 1:11
    colorBarValues(i,1) =  round(maxD - ((maxD/10)*(11-i)),2); % round(max(vq(:)) - ((max(vq(:))/10)*(11-i)),2)
    fprintf(colorBarTxt,p1Format,i,colorBarValues(i,1));
end
fclose(colorBarTxt);
cd(filePath)
% Create Heat Map as a function of Z depth
if ismember(12,outputs) == 1
    disp('7.2 Creating Shear Heat Maps for Each Frame')
    folderName = ('Through Z-Dim');
    
    delete 'HeatMaps\Shear\Through-Z Shear Stack Noise Cutoff.tif'
    delete 'HeatMaps\Shear\Through-Z Shear Stack.tif'
    
    %mkdir(strcat(filePath,'HeatMaps\'),folderName)
    for i = 1:shear.numFrames
        clear imageHeat imageHeat2 vq xq2 yq2
        [xq3,yq3] = meshgrid(0:res2*dataKey(9,1):size(Edges,2)*dataKey(9,1), 0:res2*dataKey(9,1):size(Edges,1)*dataKey(9,1));
        vq =  griddata(shear.rawX1(:),shear.rawY1(:),(shear.ltdXY(i,:)),xq3,yq3,'cubic');
        vq2 = griddata(shear.rawX1(:),shear.rawY1(:),(shear.coltdXY(i,:)),xq3,yq3,'cubic');
        xq2 = linspace(0,size(Edges,2)*dataKey(9,1),size(vq,2));
        yq2 = linspace(0,size(Edges,1)*dataKey(9,1),size(vq,1));
        
        MaximumHeatMap = imagesc(xq2,yq2,vq);
        imageHeat = MaximumHeatMap.CData;
        imageHeat(isnan(imageHeat)) = 0; 
        imageHeat = uint16(round(imageHeat * heatScale));
        imageHeatColor = ind2rgb(imageHeat,colorMap2);
        StackFile = [filePath,'\HeatMaps\Shear\','Through-Z Shear Stack.tif'];
        imwrite(imageHeatColor,StackFile,'WriteMode','append');
        
        vq2total(:,:,i) = vq2;
        MaximumHeatMap2 = imagesc(xq2,yq2,vq2);
        imageHeat2 = MaximumHeatMap2.CData;               
        imageHeat2(isnan(imageHeat2)) = 0;
        imageHeat2 = uint16(round(imageHeat2 * heatScale));
        imageHeatColor2 = ind2rgb(imageHeat2,colorMap2);
        StackFile2 = [filePath,'\HeatMaps\Shear\','Through-Z Shear Stack Noise Cutoff.tif'];
        imwrite(imageHeatColor2,StackFile2,'WriteMode','append');
        
        
        
        % %Save Max Heat Map Image
        % close all
        % maxHeatMap = figure
        % hold on
        % imshow(imageHeatColor);
        % savefile = [filePath 'HeatMaps\Through Z-Dim\HeatMapZ ' num2str(i) '.tif'];
        % if ismember(6,outputs) == 1
        %     export_fig(maxHeatMap,savefile,'-native');
        % else
        %     export_fig(maxHeatMap,savefile);
        % end
        % close
    end
    close
else
    vq2total=1;
end
end