function [imageHeatColorFinal,vqFinal] = heatmapZ(r,rDisp,planesFinal,planesGroups,imageBinary,Borders,xyScale,cutoff,colorMapZ,colormapXY,method)
filePath = strcat(cd,'\');


for i = 1:size(planesGroups,1)
    clear rVq
    maxD = 1.2; % maximum positive/negative values on Z scale bar in microns
    maxXY = 3; % maximum positive/negative values on XY scale bar in microns
    scaleD = 32768/maxD; %scalar for creating heatmapZ
    scaleXY = 65535/maxXY;
    clear planesFinal2
    planesFinal2 = planesFinal(:,planesGroups(i,1));
    if size(planesGroups,2)>1 && planesGroups(i,2)~=0
    for j = 2:nnz(planesGroups(i,:))
        planesFinal2 = cat(1,planesFinal2,planesFinal(:,planesGroups(i,j)));
    end
    end
    planesFinal2(planesFinal2==0) = [];
    
    rVq(:,1) = r(planesFinal2,1);
    rVq(:,2) = r(planesFinal2,2);
    rVq(:,3) = rDisp(planesFinal2,3);
    rVq = double(rVq);
    rVq(:,4) = sqrt(rDisp(planesFinal2,1).^2+rDisp(planesFinal2,2).^2);
    rVq(rVq(:,3)==0,:) = [];
    rVq(rVq(:,4)==0,:) = [];
    rVq(isnan(rVq(:,3)),:) = [];
    rVq(isnan(rVq(:,4)),:) = [];
    rVq(abs(rVq(:,3))<cutoff,3) = 0;
    if cutoff~=0
    rVq(abs(rVq(:,4))<.33,4) = 0;
    end
    
     
    % figure
    % scatter3(r(topMarkers,1),r(topMarkers,2),r(topMarkers,3))
    % First do Z heatmap
    res = 1;
    [xq,yq] = meshgrid(xyScale:res*xyScale:size(imageBinary,2)*xyScale, xyScale:res*xyScale:size(imageBinary,1)*xyScale);
    vq = griddata(rVq(:,1),rVq(:,2),rVq(:,3),xq,yq,'cubic');
    xq2 = linspace(0,size(imageBinary,2)*xyScale,size(vq,2));
    yq2 = linspace(0,size(imageBinary,1)*xyScale,size(vq,1));
    
    SE = strel('disk',round(5/xyScale));
    %vqEdgeFilter = isnan(vq); %imdilate(isnan(vq),SE);
    %vq = vq.*(vqEdgeFilter==0);
    if method ~=0
    vq = vq.*imageBinary.*Borders;
    end
    disp(num2str(max(max((vq)))))
    disp(num2str(min(min((vq)))))
    MaximumHeatMap = imagesc(xq2,yq2,vq);
    vqFinal(:,:,i) = single(vq);
    imageHeat = MaximumHeatMap.CData;%.*(imageBinary==0);
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
    
    savefile = [filePath strcat('HeatMaps\3D\PlanesHeatMapZ_Method_',num2str(method),'_Plane_',num2str(i),'_NoiseCutoff_',num2str(cutoff),'.tif')];
    export_fig(maxHeatMap,savefile,'-native');
    
    %Save Color Bar Values (depends on vq!)
    %Now do XY heat maps
        vq = griddata(rVq(:,1),rVq(:,2),rVq(:,4),xq,yq,'cubic');
    xq2 = linspace(0,size(imageBinary,2)*xyScale,size(vq,2));
    yq2 = linspace(0,size(imageBinary,1)*xyScale,size(vq,1));
    
    SE = strel('disk',round(5/xyScale));
    %vqEdgeFilter = isnan(vq); %imdilate(isnan(vq),SE);
    
    %vq = single(vq.*(vqEdgeFilter==0)); 
    if method ~=0
    vq = vq.*imageBinary.*Borders;
    end
    disp(num2str(max(max((vq)))))
    disp(num2str(min(min((vq)))))
    MaximumHeatMap = imagesc(xq2,yq2,vq);
    vqFinalXY(:,:,i) = single(vq);
    imageHeat = MaximumHeatMap.CData;%.*(imageBinary==0);
    imageHeat(imageHeat>0) = imageHeat(imageHeat>0)*scaleXY;
    imageHeat(isnan(imageHeat)) = 0;
    imageHeat = uint16(imageHeat);
    imageHeatColor = single(ind2rgb(imageHeat,colormapXY));
    imageHeatColorXYFinal(:,:,:,i) = imageHeatColor;
    close all
    maxHeatMap = figure;
    hold on
    imshow(imageHeatColor);
    
    savefile = [filePath strcat('HeatMaps\3D\PlanesHeatMapXY_Method_',num2str(method),'_Plane_',num2str(i),'_NoiseCutoff_',num2str(cutoff),'.tif')];
    export_fig(maxHeatMap,savefile,'-native');
    
end


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