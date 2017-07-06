function customHeatMap(book1,book2,imageBlack,dataKey,outputs,filePath)

%Input Variables to Specify
maxD = 3; %maximum value on the heat map color bar (microns)
res = .5; % number of original pixels between interpolated points (value of 1 should return a heat map with original image dimensions)
res2 = res*4;
cutoff = 2; %pixels

% Creating Heat Maps!
disp('7.1 Creating Maximum Shear Heat Maps')
close all
clear xq yq vq xq2 yq2
totalNumFrames = size(book1,2);
folderName = ('HeatMaps');
cleardir = rmdir(folderName,'s');
mkdir(filePath,folderName)
folderName = ('Single');
mkdir(strcat(filePath,'HeatMaps\'),folderName)

%Create Color Bar
colorMap2 = brewermap(65536,'*spectral');
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
savefile = [filePath '\HeatMaps\MaximumColorBar.tif'];
if ismember(6,outputs) == 1
    export_fig(colorBarSave,savefile,'-native');
else
    export_fig(colorBarSave,savefile);
end
close

%Create Heat Map of Max Displacement
[xq,yq] = meshgrid(0:res*dataKey(9,1):size(imageBlack,2)*dataKey(9,1), 0:res*dataKey(9,1):size(imageBlack,1)*dataKey(9,1));
vq = dataKey(9,1) * griddata(book2(:,1)*dataKey(9,1),book2(:,2)*dataKey(9,1),book2(:,20),xq,yq,'cubic');
vq2 = dataKey(9,1) * griddata(book2(:,1)*dataKey(9,1),book2(:,2)*dataKey(9,1),book2(:,22),xq,yq,'cubic');
disp(num2str(max(max(vq))))
xq2 = linspace(0,size(imageBlack,2)*dataKey(9,1),size(vq,2));
yq2 = linspace(0,size(imageBlack,1)*dataKey(9,1),size(vq,1));

MaximumHeatMap = imagesc(xq2,yq2,vq);
imageHeat = MaximumHeatMap.CData;
imageHeat(isnan(imageHeat)) = 0;
heatScale = (65536/maxD); %(max(max(imageHeat)))
imageHeat = uint16(round(imageHeat * heatScale));
imageHeatColor = ind2rgb(imageHeat,colorMap2);

MaximumHeatMap2 = imagesc(xq2,yq2,vq2);
imageHeat2 = MaximumHeatMap2.CData;
imageHeat(isnan(imageHeat)) = 0;
imageHeat2 = uint16(round(imageHeat2 * heatScale));
imageHeatColor2 = ind2rgb(imageHeat2,colorMap2);

%Save Max Heat Map Image
close all
maxHeatMap = figure;
hold on
imshow(imageHeatColor);

maxHeatMap2 = figure;
hold on
imshow(imageHeatColor2);

if ismember(6,outputs) == 1
    savefile = [filePath 'HeatMaps\Single\MaximumHeatMap.tif'];
    export_fig(maxHeatMap,savefile,'-native');
    savefile = [filePath 'HeatMaps\Single\MaximumHeatMapNoiseCutoff.tif'];
    export_fig(maxHeatMap2,savefile,'-native');
else
    savefile = [filePath 'HeatMaps\Single\MaximumHeatMap.tif'];
    export_fig(maxHeatMap,savefile);
    savefile = [filePath 'HeatMaps\Single\MaximumHeatMapNoiseCutoff.tif'];
    export_fig(maxHeatMap2,savefile);
end




%Save Color Bar Values (depends on vq!)
cd HeatMaps
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
    %mkdir(strcat(filePath,'HeatMaps\'),folderName)
    for i = 1:totalNumFrames
        clear imageHeat imageHeat2 vq xq2 yq2
        [xq3,yq3] = meshgrid(0:res2*dataKey(9,1):size(imageBlack,2)*dataKey(9,1), 0:res2*dataKey(9,1):size(imageBlack,1)*dataKey(9,1));
        vq = dataKey(9,1) * griddata(book2(:,1)*dataKey(9,1),book2(:,2)*dataKey(9,1),squeeze(book1(15,i,:)),xq3,yq3,'cubic');
        vq2 = dataKey(9,1) * griddata(book2(:,1)*dataKey(9,1),book2(:,2)*dataKey(9,1),squeeze(book1(19,i,:)),xq3,yq3,'cubic');
        xq2 = linspace(0,size(imageBlack,2)*dataKey(9,1),size(vq,2));
        yq2 = linspace(0,size(imageBlack,1)*dataKey(9,1),size(vq,1));
        
        MaximumHeatMap = imagesc(xq2,yq2,vq);
        imageHeat = MaximumHeatMap.CData;
        imageHeat(isnan(imageHeat)) = 0; 
        imageHeat = uint16(round(imageHeat * heatScale));
        imageHeatColor = ind2rgb(imageHeat,colorMap2);
        StackFile = [filePath,'HeatMaps\','Through-Z Shear Stack.tif'];
        imwrite(imageHeatColor,StackFile,'WriteMode','append');
        
        MaximumHeatMap2 = imagesc(xq2,yq2,vq2);
        imageHeat2 = MaximumHeatMap2.CData;               
        imageHeat2(isnan(imageHeat2)) = 0;
        imageHeat2 = uint16(round(imageHeat2 * heatScale));
        imageHeatColor2 = ind2rgb(imageHeat2,colorMap2);
        StackFile2 = [filePath,'HeatMaps\','Through-Z Shear Stack Noise Cutoff.tif'];
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
end
end