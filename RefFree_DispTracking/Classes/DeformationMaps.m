classdef DeformationMaps
    properties
        rawXY %Raw noisy XY plane deformation
        filtXY %Noise-Level Filtered XY plane deformation
        CMapXY0 %Color Map of raw data
        CMapXY %Color Map of filtered data
        
        rawZ %Raw noisy Z plane deformation
        filtZ %Noise-Level Filtered Z plane deformation
        CMapZ0 %Color Map of raw data
        CMapZ %Color Map of filtered data
    end
    
    methods
        %% Spatial Maps of Deformation Magnitude
        function obj = DeformationMaps(r,m,raw,plane,image)            
             pSize = raw.dataKey(9,1);
            for i = 1:size(plane.groups,1)
                %  Prepare indexed list of features in current plane
                %Some 'planes' are separated, use plane.groups variable to concatenate
                %them
                clear AllFeatures
                AllFeatures = plane.final(:,plane.groups(i,1));
                if size(plane.groups,2)>1 && plane.groups(i,2)~=0
                    for j = 2:nnz(plane.groups(i,:))
                        AllFeatures = cat(1,AllFeatures,plane.final(:,plane.groups(i,j)));
                    end
                end
                AllFeatures(AllFeatures==0) = []; %Clean list of features
                
                % Write XYZ displacement data of current plane to a list
                clear rVq
                rVq(:,1) = r.X(AllFeatures);
                rVq(:,2) = r.Y(AllFeatures);
                rVq(:,3) = r.dZ(AllFeatures);
                rVq(:,4) = r.dS(AllFeatures);
                rVq = double(rVq);
                
                % Remove incomplete datapoints
                rVq(rVq(:,4)==0,:) = []; %eliminates features where xy information is not usable
                rVq(rVq(:,4)>7,:) = []; %eliminates features where xy information is not usable
                rVq(isnan(rVq(:,3)),:) = [];
                rVq(isnan(rVq(:,4)),:) = [];
                
                % Check to see if any points remain, if not:
                if size(rVq,1) == 0
                    rVq(1,1:4) = 0;
                end                
                res = 2.12/0.1625;
                [xq,yq] = meshgrid(pSize:res*pSize:size(image.ADil,2)*pSize, pSize:res*pSize:size(image.ADil,1)*pSize);
                
                % Map of Z displacements
                vq = griddata(rVq(:,1),r.s(1,2)-rVq(:,2),rVq(:,3),xq,yq,'cubic');                                
                obj.rawZ(:,:,i) = single(vq);
                vq = vq.*flipud(imresize(image.imgNBds.*image.Borders,size(vq))>0);
                vq(abs(vq)<m.ZnoiseCO) = 0;
                obj.filtZ(:,:,i) = single(vq);
                
                % Map of Shear displacement
                vq = griddata(rVq(:,1),r.s(1,2)-rVq(:,2),rVq(:,4),xq,yq,'cubic');                                
                obj.rawXY(:,:,i) = single(vq);
                vq = vq.*flipud(imresize(image.imgNBds.*image.Borders,size(vq))>0);
                vq(abs(vq)<m.SnoiseCO) = 0;
                obj.filtXY(:,:,i) = single(vq);
            end
        end
        %% Colorized Maps of Deformation
        function obj = ZDeformationColorMap(obj,image,m,raw,maxD,cmapType)
            pSize = raw.dataKey(9,1);
            cmap = single(brewermap(65536,cmapType));
            filePath = strcat(cd,'\');
            scaleD = 32768/maxD;
            xq2 = linspace(0,size(image.ADil,2)*pSize,size(obj.rawZ,2));
            yq2 = linspace(0,size(image.ADil,1)*pSize,size(obj.rawZ,1));
            for i = 1:size(obj.rawZ,3)
                MaximumHeatMap = imagesc(xq2,yq2,obj.rawZ(:,:,i));
                imageHeat = MaximumHeatMap.CData;%.*(image.ADil==0);
                imageHeat = imresize(imageHeat,size(image.ADil),'nearest');              
                imageHeat(imageHeat>0) = 32768+(abs(imageHeat(imageHeat>0))*scaleD);
                imageHeat(imageHeat<0) = 32768 - (abs(imageHeat(imageHeat<0))*scaleD);
                imageHeat(imageHeat==0) = 32768;
                imageHeat(isnan(imageHeat)) = 32768;
                imageHeat = uint16(imageHeat);
                imageHeatColor = single(ind2rgb(imageHeat,cmap));
                obj.CMapZ0(:,:,:,i) = imageHeatColor;
                close all
                maxHeatMap = figure;
                hold on
                imshow(imageHeatColor);
                savefile = [filePath strcat('HeatMaps\3D\PlanesHeatMapZ_Method_',num2str(0),'_Plane_',num2str(i),'_NoiseCutoff_',num2str(0),'.tif')];
                export_fig(maxHeatMap,savefile,'-native');

                MaximumHeatMap = imagesc(xq2,yq2,obj.filtZ(:,:,i));
                imageHeat = MaximumHeatMap.CData;%.*(image.ADil==0);            
                imageHeat = imresize(imageHeat,size(image.ADil),'bicubic');                
                imageHeat(imageHeat>0) = 32768+(abs(imageHeat(imageHeat>0))*scaleD);
                imageHeat(imageHeat<0) = 32768 - (abs(imageHeat(imageHeat<0))*scaleD);
                imageHeat(imageHeat==0) = 32768;
                imageHeat(isnan(imageHeat)) = 32768;
                imageHeat = uint16(imageHeat);
                imageHeatColor = single(ind2rgb(imageHeat,cmap));
                obj.CMapZ(:,:,:,i) = imageHeatColor;
                close all
                maxHeatMap = figure;
                hold on
                imshow(imageHeatColor);
                savefile = [filePath strcat('HeatMaps\3D\PlanesHeatMapZ_Method_',num2str(3),'_Plane_',num2str(i),'_NoiseCutoff_',num2str(m.ZnoiseCO),'.tif')];
                export_fig(maxHeatMap,savefile,'-native');
            end
            
            mkdir('HeatMaps\3D\ColorBar')            
            colorBar1 = single(zeros(500,25));
            range = uint16(round(linspace(65536,1,500)'));
            for i = 1:25
                colorBar1(1:500,i) = range;
            end
            colorBar2 = ind2rgb(colorBar1,cmap);
            for i = 1:10
                colorBar2((i*50)-3:(i*50),13:25,:) = 0;
            end
            colorBar2(1:3,13:25,:) = 0;
            
            %Save Color Bar Image
            close all
            colorBarSave = figure;
            hold on
            imshow(colorBar2);
            filePath=cd;
            savefile = [filePath '\HeatMaps\3D\ColorBar\ColorBarZ.tif'];
            export_fig(colorBarSave,savefile,'-native');
        end
        %%
        function obj = XYDeformationColorMap(obj,image,m,raw,maxD,cmapType)
            pSize = raw.dataKey(9,1);
            cmap = single(brewermap(65536,cmapType));
            filePath = strcat(cd,'\');
            scaleXY = 65535/maxD;
            xq2 = linspace(0,size(image.ADil,2)*pSize,size(obj.rawXY,2));
            yq2 = linspace(0,size(image.ADil,1)*pSize,size(obj.rawXY,1));
            for i = 1:size(obj.rawXY,3)
                %Raw Map                
                MaximumHeatMap = imagesc(xq2,yq2,obj.rawXY(:,:,i));
                imageHeat = MaximumHeatMap.CData;%.*(image.ADil==0);
                imageHeat = imresize(imageHeat,size(image.ADil),'nearest');
                imageHeat(imageHeat>0) = imageHeat(imageHeat>0)*scaleXY;
                imageHeat(isnan(imageHeat)) = 0;
                imageHeat = uint16(imageHeat);
                imageHeatColor = single(ind2rgb(imageHeat,cmap));
                obj.CMapXY0(:,:,:,i) = imageHeatColor;
                close all
                maxHeatMap = figure;
                hold on
                imshow(imageHeatColor);
                savefile = [filePath strcat('HeatMaps\3D\PlanesHeatMapXY_Method_',num2str(0),'_Plane_',num2str(i),'_NoiseCutoff_',num2str(0),'.tif')];
                export_fig(maxHeatMap,savefile,'-native');
                
                %Filtered Map
                MaximumHeatMap = imagesc(xq2,yq2,obj.filtXY(:,:,i));
                imageHeat = MaximumHeatMap.CData;%.*(image.ADil==0);
                imageHeat = imresize(imageHeat,size(image.ADil),'bicubic');               
                imageHeat(imageHeat>0) = imageHeat(imageHeat>0)*scaleXY;
                imageHeat(isnan(imageHeat)) = 0;
                imageHeat = uint16(imageHeat);
                imageHeatColor = single(ind2rgb(imageHeat,cmap));
                obj.CMapXY(:,:,:,i) = imageHeatColor;
                close all
                maxHeatMap = figure;
                hold on
                imshow(imageHeatColor);
                
                savefile = [filePath strcat('HeatMaps\3D\PlanesHeatMapXY_Method_',num2str(3),'_Plane_',num2str(i),'_NoiseCutoff_',num2str(m.SnoiseCO),'.tif')];
                export_fig(maxHeatMap,savefile,'-native');
            end
            
            %Make Color Bar
            mkdir('HeatMaps\3D\ColorBar')            
            colorBar1 = single(zeros(500,25));
            range = uint16(round(linspace(65536,1,500)'));
            for i = 1:25
                colorBar1(1:500,i) = range;
            end
            colorBar2 = single(ind2rgb(colorBar1,cmap));
            for i = 1:10
                colorBar2((i*50)-3:(i*50),13:25,:) = 0;
            end
            colorBar2(1:3,13:25,:) = 0;
            
            %Save Color Bar Image
            close all
            colorBarSave = figure;
            hold on
            imshow(colorBar2);
            filePath=cd;
            savefile = [filePath '\HeatMaps\3D\ColorBar\ColorBarShear.tif'];
            export_fig(colorBarSave,savefile,'-native');
                        
                
        end
        %%
        % Calculate Useful parameters
        function [obj]  = vqStats(obj,plane)
            vq2pos = vq2(:,:,:).*(vq2(:,:,:)>0);
            vq2pos(isnan(vq2pos)) = 0;
            vq2neg = vq2(:,:,:).*(vq2(:,:,:)<0);
            vq2neg(isnan(vq2neg)) = 0;
            for i = 1:size(plane.groups,1)
                rDisp2PosTotal(i,1) = sum(sum(vq2pos(:,:,i)));
                rDisp2PosMax(i,1) = max(max(vq2pos(:,:,i)));
                rDisp2NegTotal(i,1) = sum(sum(vq2neg(:,:,i)));
                rDisp2NegMax(i,1) = max(max(abs(vq2neg(:,:,i)))) * -1;
            end
        end
    end
end