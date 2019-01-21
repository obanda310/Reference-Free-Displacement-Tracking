function [imageHeatColorFinal,vqFinal] = heatmapZ(r,rDisp,planesFinal,planesGroups,imageBinary,Borders,xyScale,cutoff,colorMapZ,colormapXY,method)

% This function is designed to create a heat map of both shear and normal
% displacements. There are several versions of these heatmaps which are
% generated. A heat map is created for every registered plane of features
% in the planesGroups variable. For each plane, 3 versions of each heatmap
% are generated: Full data, data within imageBinary, and data outside of
% imageBinary. Each has had a purpose so far.

%% Set global limits and constants

filePath = strcat(cd,'\');
maxD = 1.2; % maximum positive/negative values on Z scale bar in microns
maxXY = 3; % maximum positive/negative values on XY scale bar in microns
scaleD = 32768/maxD; %scalar for creating heatmapZ
scaleXY = 65535/maxXY;

for i = 1:size(planesGroups,1)
    %%  Prepare indexed list of features in current plane
    %Some 'planes' are separated, use planesGroups variable to concatenate
    %them
    clear AllFeatures
    AllFeatures = planesFinal(:,planesGroups(i,1));
    if size(planesGroups,2)>1 && planesGroups(i,2)~=0
        for j = 2:nnz(planesGroups(i,:))
            AllFeatures = cat(1,AllFeatures,planesFinal(:,planesGroups(i,j)));
        end
    end
    AllFeatures(AllFeatures==0) = []; %Clean list of features
    
    %% Write XYZ displacement data of current plane to a list
    clear rVq
    rVq(:,1) = r(AllFeatures,1);
    rVq(:,2) = r(AllFeatures,2);
    rVq(:,3) = rDisp(AllFeatures,3);
    rVq(:,4) = sqrt(rDisp(AllFeatures,1).^2+rDisp(AllFeatures,2).^2);
    rVq = double(rVq);
    
    %% Remove incomplete datapoints
    rVq(rVq(:,4)==0,:) = []; %eliminates all features where xy information is not usable
    rVq(isnan(rVq(:,3)),:) = [];
    rVq(isnan(rVq(:,4)),:) = [];
    rVq(abs(rVq(:,3))<cutoff,3) = 0;
    if cutoff~=0
        rVq(abs(rVq(:,4))<.33,4) = 0;
    end
    % Check to see if any points remain, if not:
    if size(rVq,1) == 0
        rVq(1,1:4) = 0;
    end
    
    averageZ = round(mean(r(AllFeatures,3)),2);
    %% Heatmap of Z displacements
    
    res = 2.12/0.1625;
    [xq,yq] = meshgrid(xyScale:res*xyScale:size(imageBinary,2)*xyScale, xyScale:res*xyScale:size(imageBinary,1)*xyScale);
    vq = griddata(rVq(:,1),rVq(:,2),rVq(:,3),xq,yq,'cubic');
    xq2 = linspace(0,size(imageBinary,2)*xyScale,size(vq,2));
    yq2 = linspace(0,size(imageBinary,1)*xyScale,size(vq,1));
    
    if cutoff ~= 0
        vq = vq.*(imresize(imageBinary.*Borders,size(vq))>0);
    end
    %     disp(num2str(max(max((vq)))))
    %     disp(num2str(min(min((vq)))))
    MaximumHeatMap = imagesc(xq2,yq2,vq);
    vqFinal(:,:,i) = single(vq);
    imageHeat = MaximumHeatMap.CData;%.*(imageBinary==0);
    if method == 3
        imageHeat = imresize(imageHeat,size(imageBinary),'bicubic');
    else
        imageHeat = imresize(imageHeat,size(imageBinary),'nearest');
    end
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
    
    savefile = [filePath strcat('HeatMaps\3D\PlanesHeatMapZ_Method_',num2str(method),'_Plane_',num2str(averageZ),'_NoiseCutoff_',num2str(cutoff),'.tif')];
    export_fig(maxHeatMap,savefile,'-native');
    
    %% Heat map of Shear displacement
    vq = griddata(rVq(:,1),rVq(:,2),rVq(:,4),xq,yq,'cubic');
    xq2 = linspace(0,size(imageBinary,2)*xyScale,size(vq,2));
    yq2 = linspace(0,size(imageBinary,1)*xyScale,size(vq,1));
    
    if method ~=0
        vq = vq.*(imresize(imageBinary.*Borders,size(vq))>0);
    end
    %     disp(num2str(max(max((vq)))))
    %     disp(num2str(min(min((vq)))))
    MaximumHeatMap = imagesc(xq2,yq2,vq);
    vqFinalXY(:,:,i) = single(vq);
    imageHeat = MaximumHeatMap.CData;%.*(imageBinary==0);
    if method == 3
        imageHeat = imresize(imageHeat,size(imageBinary),'bicubic');
    else
        imageHeat = imresize(imageHeat,size(imageBinary),'nearest');
    end
    imageHeat(imageHeat>0) = imageHeat(imageHeat>0)*scaleXY;
    imageHeat(isnan(imageHeat)) = 0;
    imageHeat = uint16(imageHeat);
    imageHeatColor = single(ind2rgb(imageHeat,colormapXY));
    imageHeatColorXYFinal(:,:,:,i) = imageHeatColor;
    close all
    maxHeatMap = figure;
    hold on
    imshow(imageHeatColor);
    
    savefile = [filePath strcat('HeatMaps\3D\PlanesHeatMapXY_Method_',num2str(method),'_Plane_',num2str(averageZ),'_NoiseCutoff_',num2str(cutoff),'.tif')];
    export_fig(maxHeatMap,savefile,'-native');
    
end

%% Create a text document with color bar values
cd HeatMaps\3D\ColorBar
colorBarTxt = fopen('Heat Map Color Bar Values Z.txt','wt');
p1Format = 'Tick number %1.0f is %.2f \n';
for i = 1:11
    colorBarValues(i,1) =  round(maxD - ((maxD/5)*(11-i)),2); % round(max(vq(:)) - ((max(vq(:))/10)*(11-i)),2)
    fprintf(colorBarTxt,p1Format,i,colorBarValues(i,1));
end
fclose(colorBarTxt);
cd(filePath)

end