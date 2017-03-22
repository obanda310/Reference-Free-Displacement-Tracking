function customHeatMap(book1,book2,imageBlack,dataKey,outputs,filePath)
% Creating Heat Maps!
disp('Creating Heat Maps')
close all
clear xq yq vq xq2 yq2
totalNumFrames = size(book1,2);
folderName = ('HeatMaps');
cleardir = rmdir(folderName,'s');
mkdir(filePath,folderName)
folderName = ('Single');
mkdir(strcat(filePath,'HeatMaps\'),folderName)

%Create Color Bar
colorMap2 = brewermap(65536,'*Spectral');
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
colorBarSave = figure
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
[xq,yq] = meshgrid(0:2*dataKey(9,1):size(imageBlack,2)*dataKey(9,1), 0:2*dataKey(9,1):size(imageBlack,1)*dataKey(9,1));
vq = dataKey(9,1) * griddata(book2(:,1)*dataKey(9,1),book2(:,2)*dataKey(9,1),book2(:,13),xq,yq);
xq2 = linspace(0,size(imageBlack,2)*dataKey(9,1),size(vq,2));
yq2 = linspace(0,size(imageBlack,1)*dataKey(9,1),size(vq,1));
MaximumHeatMap = imagesc(xq2,yq2,vq);
imageHeat = MaximumHeatMap.CData;
imageHeat(isnan(imageHeat)) = 0;
heatScale = (65536/2.25); %(max(max(imageHeat)))
imageHeat = uint16(round(imageHeat * heatScale));
imageHeatColor = ind2rgb(imageHeat,colorMap2);

%Save Max Heat Map Image
close all
maxHeatMap = figure
hold on
imshow(imageHeatColor);
savefile = [filePath 'HeatMaps\Single\MaximumHeatMap.tif'];
if ismember(6,outputs) == 1
    export_fig(maxHeatMap,savefile,'-native');
else
    export_fig(maxHeatMap,savefile);
end
close

%Save Color Bar Values (depends on vq!)
cd HeatMaps
colorBarTxt = fopen('Heat Map Color Bar Values.txt','wt');
p1Format = 'Tick number %1.0f is %.2f \n';
for i = 1:11
    colorBarValues(i,1) =  round(2.25 - ((2.25/10)*(11-i)),2); % round(max(vq(:)) - ((max(vq(:))/10)*(11-i)),2)
    fprintf(colorBarTxt,p1Format,i,colorBarValues(i,1));
end
fclose(colorBarTxt);
cd(filePath)
% Create Heat Map as a function of Z depth
if ismember(12,outputs) == 1
    folderName = ('Through Z-Dim');
    %mkdir(strcat(filePath,'HeatMaps\'),folderName)
    for i = 1:totalNumFrames
        clear imageHeat vq xq2 yq2
        
        vq = dataKey(9,1) * griddata(book2(:,1)*dataKey(9,1),book2(:,2)*dataKey(9,1),squeeze(book1(10,i,:)),xq,yq);
        xq2 = linspace(0,size(imageBlack,2)*dataKey(9,1),size(vq,2));
        yq2 = linspace(0,size(imageBlack,1)*dataKey(9,1),size(vq,1));
        MaximumHeatMap = imagesc(xq2,yq2,vq);
        imageHeat = MaximumHeatMap.CData;
        imageHeat(isnan(imageHeat)) = 0;
        imageHeat = uint16(round(imageHeat * heatScale));
        imageHeatColor = ind2rgb(imageHeat,colorMap2);
        StackFile = [filePath,'HeatMaps\','NormalStack.tif'];
        imwrite(imageHeatColor,StackFile,'WriteMode','append');
        
        
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
end
% Create a GAUSS filtered version of a Heat Map as a function of Z depth
if ismember(13,outputs) == 1
    folderName = ('Through Z-Dim Gauss');
    %mkdir(strcat(filePath,'HeatMaps\'),folderName)
    for i = 1:totalNumFrames
        clear imageHeat vq xq2 yq2
        
        vq = dataKey(9,1) * griddata(book2(:,1)*dataKey(9,1),book2(:,2)*dataKey(9,1),squeeze(book1(10,i,:)),xq,yq);
        vqGauss = imgaussfilt(vq,round((1/dataKey(9,1))*(size(vq,2)/size(imageBlack,2))));
        xq2 = linspace(0,size(imageBlack,2)*dataKey(9,1),size(vq,2));
        yq2 = linspace(0,size(imageBlack,1)*dataKey(9,1),size(vq,1));
        MaximumHeatMap = imagesc(xq2,yq2,vqGauss);
        imageHeat = MaximumHeatMap.CData;
        imageHeat(isnan(imageHeat)) = 0;
        imageHeat = uint16(round(imageHeat * heatScale));
        imageHeatColor = ind2rgb(imageHeat,colorMap2);
        GaussStackFile = [filePath,'HeatMaps\','GaussStack.tif'];
        imwrite(imageHeatColor,GaussStackFile,'WriteMode','append');
        
        % %Save Max Heat Map Image
        % close all
        % maxHeatMap = figure
        % hold on
        % imshow(imageHeatColor);
        % savefile = [filePath 'HeatMaps\Through Z-Dim Gauss\HeatMapZ ' num2str(i) '.tif'];
        % if ismember(6,outputs) == 1
        %     export_fig(maxHeatMap,savefile,'-native');
        % else
        %     export_fig(maxHeatMap,savefile);
        % end
        % close
    end
end
close
end