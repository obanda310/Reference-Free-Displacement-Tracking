classdef RawData3D
    properties
        s
        r
        X
        Y
        Z
        rX %reference coordinates
        rY %reference coordinates
        rZ %reference coordinates
        dX %displacement from reference
        dY %displacement from reference
        dZ %displacement from reference
        dS %total displacement in XY plane
        l
        Mass
        MaxI
        rg
        PctAbove
        row %8
        rVa % row end 1
        rVb % row end 2
        col %9
        colSize
        colUp %the index of the next object upwards in the column (if it exists)
        colDown
        rowRight
        rowLeft
        plane
        planeGroup
        vplane
        colCheck %logical value describing whether the next object in this objects column exists in the next plane (i.e. the next object isn't missing.)
        XSC %10 Shift correction X
        YSC %11 Shift Correction Y
        D %markers in deformation zone
        Dl %logical version of D
        ND %non-deformed markers
        NDl %logical version of ND
        State %state of maturation of the instance of rawdata3d (increases as the main script runs)
        
    end
    methods
        %%
        function obj = RawData3D(res,raw,r)
            if nargin == 2
                obj.s(1,1) = size(res,1); obj.s(2,1) = size(res,2); obj.s(3,1) = size(res,3);
                obj.s(1:2,2) = obj.s(1:2,1)*raw.dataKey(9,1);
                obj.s(4,2) = obj.s(3,1)*raw.dataKey(10,1);
                d = round(0.8/raw.dataKey(9,1));
                dz = round(2/raw.dataKey(10,1));
                obj.r=...
                    feature3dMB(res, d , [d d dz], [obj.s(1) obj.s(2) obj.s(3)],[1 1 1],round(.5/raw.dataKey(9,1)),0,.8,.9); %
                obj.r(:,1:2) = obj.r(:,1:2)*raw.dataKey(9,1);
                
                % remove detections near edges
                borderwidth = 1;% should be slightly larger than half the XY radius of an object
                tempr = obj.r;
                tempr = floor(tempr);
                removeidx = tempr(:,1)<=borderwidth | tempr(:,2)<=borderwidth | obj.s(1,2)-tempr(:,1)<=borderwidth | obj.s(2,2)-tempr(:,2)<=borderwidth;
                obj.r(removeidx,:) = [];
                removeidx = obj.r(:,4)<2000;
                obj.r(removeidx,:) = [];
                %adjust X and Y data to match image orientation when
                %plotted (invert X, then switch X and Y)
                obj.r(:,1) = obj.s(1,2)-obj.r(:,1);
                temp = obj.r(:,1);
                obj.r(:,1) = obj.r(:,2);
                obj.r(:,2) = temp;
                
                obj.r(:,3) = obj.r(:,3)*raw.dataKey(10,1);
                obj.l = size(obj.r,1);
            elseif nargin == 3
                obj.s(1,1) = size(res,1); obj.s(2,1) = size(res,2); obj.s(3,1) = size(res,3);
                obj.s(1:2,2) = obj.s(1:2,1)*raw.dataKey(9,1);
                obj.r = r;
                obj.l = size(r,1);
            end
        end
        %%
        function obj = TranscribeR(obj)
            obj.X = obj.r(:,1);
            obj.Y = obj.r(:,2);
            obj.Z = obj.r(:,3);
            if size(obj.r,2)>3
                obj.Mass = obj.r(:,4);
                obj.rg = obj.r(:,5);
                obj.MaxI = obj.r(:,6);
                obj.PctAbove = obj.r(:,7);
            end
        end
        %%
        function viewDetections(obj,raw,offset,image)
            figure
            if nargin == 2
                scatter3(obj.X,obj.Y,obj.Z)
            elseif nargin == 4
                imshow(image)
                hold on
                scatter3((obj.X-offset(1,1))/raw.dataKey(9,1),((obj.Y)-offset(1,2))/raw.dataKey(9,1),obj.Z)
                hold off
            end
                        bcolor = 'black';
            fcolor = 'white';
            AxisFontSize = 12;
            LegendFontSize = 14;
            xt = 'X \mum';% input('enter the xaxis label','s');
            yt = 'Y \mum'; %input('enter the yaxis label','s');
            zt = 'Z \mum';
            label{1} = xlabel(xt);
            label{2} = ylabel(yt);
            label{3} = zlabel(zt);
            set(gca,'YMinorTick','on','color',bcolor)
            ytickformat('%.1f')
            le{1} = 'plane 1'; %input('enter the legend','s');
            le{2} = 'plane 2';
            ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,1,[0 0])
            %errorbar(meanDisplacements(1,1:3),meanDisplacements(2,1:3),'.','color',[0 0 0],'MarkerSize',1)
            
            legend off
            hold off
            view(-41,33)
        end
        %%
        function obj = regionCheck(obj,imageADil,raw)
            obj.ND = zeros(1,1);
            obj.D = zeros(1,1);
            obj.Dl = zeros(size(obj.X));
            obj.NDl = zeros(size(obj.X));
            for i = 1:obj.l
                %if it is under the cell
                if (obj.s(1,2)-round(obj.Y(i))) <=0 || round(obj.X(i)/raw.dataKey(9,1)) <= 0 || (obj.s(1,2)-round(obj.Y(i)))/raw.dataKey(9,1) >= obj.s(1,1) || round(obj.X(i)/raw.dataKey(9,1)) >= obj.s(2,1)
                    obj.ND = cat(1,obj.ND,i);
                    obj.NDl(i,1) = 1;
                elseif imageADil(obj.s(1,1)-(round(obj.Y(i)/raw.dataKey(9,1))),round(obj.X(i)/raw.dataKey(9,1)))~=0
                    obj.ND = cat(1,obj.ND,i);
                    obj.NDl(i,1) = 1;
                else
                    obj.D = cat(1,obj.D,i);
                end
            end
            obj.ND(1,:) = [];
            obj.D(1,:) = [];
            obj.Dl = obj.NDl==0;
        end
        %% Augment Non-Deformed Marker list with row ends
        function obj = AugmentND(obj,rows)
         % Grab objects on the edge
            tic
            tidx2=[];
            tidx = 1:1:obj.l;
            [pIdx,~] = dsearchn(obj.r(tidx,1:3),obj.r(tidx,1:3)+(2*rows.V(1:3,1))');
            for i = 1:nnz(tidx)
                if tidx(pIdx(i)) == tidx(i)
                    tidx2 = cat(1,tidx2,tidx(i));
                end     
            end
            
            [pIdx,~] = dsearchn(obj.r(tidx,1:3),obj.r(tidx,1:3)-(2*rows.V(1:3,1))');
            for i = 1:nnz(tidx)
                if tidx(pIdx(i)) == tidx(i)
                    tidx2 = cat(1,tidx2,tidx(i));
                end     
            end
                                    
            for i = 1:nnz(tidx2)
                obj.NDl(tidx2(i)) = 1;
            end
            
            obj.NDl = obj.NDl>0;
            obj.ND = find(obj.NDl);
            obj.Dl = obj.NDl==0;
            obj.D = find(obj.Dl);
            
        end
        %% Procedurally Augment Non-Deformed Marker list with objects which are probably not deformed
        function obj = PAugmentND(obj,rows)
            
            
        end
        %%
        function testRegionCheck(obj,image)
            figure
            imshow(image.ADil,[])
            
            figure
            hold on
            for i = obj.ND
                scatter3(obj.X(i),obj.Y(i),obj.Z(i))
            end
            
            clear i
            figure
            hold on
            for i = [1:size(obj.X,1)]'
                scatter3(obj.X(i),obj.Y(i),obj.Z(i))
            end
        end
        %%
        function obj = updatePlane(obj,plane)
            obj.plane = zeros(size(obj.X));
            for i = 1:size(plane.final,2)
                for j = 1:size(plane.final,1)
                    if plane.final(j,i) > 0
                        obj.plane(plane.final(j,i)) = i;
                    end
                end
            end
            obj.planeGroup = zeros(size(obj.X));
            for i = 1:size(obj.r,1)
                [tempx,tempy] = find(plane.groups == obj.plane(i));
                obj.planeGroup(i,1) = tempx;
            end
        end
        %%
        function obj = updateColumn(obj)
            if obj.State>1
            rFlat(:,1) = obj.rX;
            rFlat(:,2) = obj.rY;            
            else
            rFlat(:,1) = obj.X;
            rFlat(:,2) = obj.Y;            
            end
            rFlat(:,6) = max(obj.planeGroup)-obj.planeGroup+1;
            
            % Here we will link objects in the Z direction
            disp('Linking 3D Detections Between Planes.')
            % Set a maximum linking distance in microns that any object can still be
            % considered part of a pillar. Smaller values will speed up code.
            maxLinkDistance = .9; % !!**value should be smaller than array spacing**!!
            maxLD = maxLinkDistance;% /raw.dataKey(9,1); <--if using pixels, use this conversion
            disp(['Maximum allowable XY Translation (Microns): ',num2str(maxLinkDistance)])
            % Set a maximum number of frames to look for a linked object before giving
            % up (maxJumpDistance)
            maxJD = 1;
            disp(['Max Jump Distance (Number of Possible Empty Planes between detections): ',num2str(maxJD)])
            rFlat(:,3) = 1:size(rFlat,1); %record pillar indices for posterity
            [rFlat2,order] = sortrows(rFlat,6); %column six needs to be in increasing order for trackmem to work
            tidx = rFlat2(:,1)==0 | rFlat2(:,2)==0;
            rFlat2(tidx,:) = [];
            % Linking using trackmem function (kilfoil code)
            disp('Linking Pillars')
            [lub] = trackmem(rFlat2,maxLD,2,0,maxJD);
            % remove bad links (will only work if non-deformed columns can be linked accurately)
            
            % **First, establish the average Non-Deformed column**
            % keep only non-deformed objects
            lub2 = lub;
            lub2 = sortrows(lub2,3);
            for i = 1:size(lub2,1)
                if ismember(lub2(i,3),obj.ND)==0
                    lub2(i,:) = NaN;
                end
            end
            lub2(isnan(lub2(:,1)),:) = [];
            lub2 = sortrows(lub2,[7 6]);
            
            % remove columns with fewer than max planes            
            for i = 1:max(lub2(:,7))
                lub2(lub2(:,7)==i,4) = size(lub2(lub2(:,7)==i,1),1);
                if size(lub2(lub2(:,7)==i,1),1) == 1
                    lub2(lub2(:,7)==i,7) = 0;
                else
                    lub2(lub2(:,7)==i,4) = size(lub2(lub2(:,7)==i,1),1);
                end
            end
            lub2(lub2(:,4)<max(lub2(:,4)),:) = [];
            
            for i = 1:size(lub2,1)-1
                if lub2(i,7) == lub2 (i+1,7)
                    lub2(i+1,4) = lub2(i+1,1) -lub2(i,1);
                    lub2(i+1,5) = lub2(i+1,2) -lub2(i,2);
                else
                    lub2(i+1,4) = 0;
                    lub2(i+1,5) = 0;
                end
            end
            
            %store  the average shift in successive columns to be used as a
            %correction later
            adj = zeros(max(lub2(:,6)),2);
            for i = 1:max(lub2(:,6))
                idx = lub2(:,6) == i;
                adj(i,1) = mean(lub2(idx,4));
                adj(i,2) = mean(lub2(idx,5));
            end
            
            % Then calculate the differences between successive objects in
            % column, adjusting for the average column
            for i = 1:size(lub,1)-1
                if lub(i,7) == lub (i+1,7)
                    tplane = lub(i+1,6);
                    lub(i+1,4) = (lub(i+1,1)-lub(i,1))- adj(tplane,1);
                    lub(i+1,5) = (lub(i+1,2)-lub(i,2))- adj(tplane,2);
                else
                    lub(i+1,4) = 0;
                    lub(i+1,5) = 0;
                end
            end
            
            
            for i = 2:size(lub,1)
                if (lub(i,7) == lub (i-1,7)) && (lub(i-1,4) > 0.2 )
                    if lub(i,4) < lub(i-1,4)
                        lub(i,7) = max(lub(:,7))+1;
                    end
                elseif (lub(i,7) == lub (i-1,7)) && (lub(i-1,4) < -0.2 )
                    if lub(i,4) > lub(i-1,4)
                        lub(i,7) = max(lub(:,7))+1;
                    end
                elseif (lub(i,7) == lub (i-1,7)) &&  lub(i-1,5) > 0.2
                    if lub(i,5) < lub(i-1,5)
                        lub(i,7) = max(lub(:,7))+1;
                    end
                elseif (lub(i,7) == lub (i-1,7)) &&  lub(i-1,5) < -0.2
                    if lub(i,5) > lub(i-1,5)
                        lub(i,7) = max(lub(:,7))+1;
                    end
                end
            end
            lub = sortrows(lub,[7,6]);
            
            
            % Store pillar size
            for i = 1:max(lub(:,7))
                lub(lub(:,7)==i,4) = size(lub(lub(:,7)==i,1),1);
                if size(lub(lub(:,7)==i,1),1) == 1 && (lub(lub(:,7)==i,6)~=min(lub(:,6)))
                    lub(lub(:,7)==i,7) = 0;                
                else
                    lub(lub(:,7)==i,4) = size(lub(lub(:,7)==i,1),1);
                end
            end
            lub = sortrows(lub,3);
            % Update rawData variable with column info
            obj.col = zeros(size(obj.X));
            obj.colSize = zeros(size(obj.X));
            obj.col(lub(:,3),1) = lub(:,7);
            obj.colSize(lub(:,3),1) = lub(:,4);
            
        end
        
        %%
        function obj = updateVPlane(obj,VPlanes,preview)
            obj.vplane = zeros(size(obj.X));
            for i = 1:size(VPlanes.rows,1)
                for j = 1:size(VPlanes.rows,2)
                    if VPlanes.rows(i,j)>0
                        tIdx = obj.row == VPlanes.rows(i,j);
                        obj.vplane(tIdx,1) = i;
                    end
                end
            end
            for i = 1:size(VPlanes.cols,1)
                for j = 1:size(VPlanes.cols,2)
                    if VPlanes.cols(i,j)>0
                        tIdx = obj.col == VPlanes.cols(i,j);
                        obj.vplane(tIdx,1) = i;
                    end
                end
                
            end
            
            if preview >0
                figure
                hold on
                for i = 1:max(obj.vplane)
                    tIdx = obj.vplane == i;
                    scatter3(obj.X(tIdx),obj.Y(tIdx),obj.Z(tIdx))
                end
            end
            
        end
        %% Update reference values based on rowfits
        function obj = updateReference(obj,m3)
            obj.rX = m3.ref(:,1);
            obj.rY = m3.ref(:,2);
            obj.rZ = m3.ref(:,3);
            
            obj.XSC = [];
            obj.YSC = [];
            
            
            obj.XSC = cat(1,m3.refSC(:,1),obj.XSC);
            obj.YSC = cat(1,m3.refSC(:,2),obj.YSC);
            
            
            obj.rVa = zeros(obj.l,3);
            obj.rVb = zeros(obj.l,3);
            for i = 1:obj.l
                if obj.row(i) >0
                    obj.rVa(i,1:3) = m3.rowFits(obj.row(i),1:3,1);
                    obj.rVb(i,1:3) = m3.rowFits(obj.row(i),1:3,2);
                end
            end
            obj.dX = zeros(size(obj.rX));
            obj.dY = zeros(size(obj.rY));
            obj.dZ = zeros(size(obj.rZ));
            obj.dS = zeros(size(obj.X));
            clear tidx
            tidx = obj.rX>0;
            obj.dX(tidx) = obj.X(tidx)-obj.rX(tidx);
            obj.dY(tidx) = obj.Y(tidx)-obj.rY(tidx);
            obj.dZ(tidx) = obj.Z(tidx)-obj.rZ(tidx);
            obj.dS(tidx) = sqrt(obj.dX(tidx).^2+obj.dY(tidx).^2);
        end
        %% Reset certain objects' data
        function obj = resetData(obj,idx)
            for i =1:size(idx,1)
                obj.colUp(obj.colUp == idx(i))=0;
                obj.colDown(obj.colDown == idx(i))=0;
                obj.rowLeft(obj.rowLeft== idx(i))=0;
                obj.rowRight(obj.rowRight== idx(i))=0;
            end
            obj.rX(idx) = 0;
            obj.rY(idx) = 0;
            obj.rZ(idx) = 0;
            obj.dX(idx) = 0;
            obj.dY(idx) = 0;
            obj.dZ(idx) = 0;
            obj.dS(idx) = 0;
            obj.row(idx) = 0;
            obj.rVa(idx) = 0;
            obj.rVb(idx) = 0;
            obj.col(idx) = 0;
            obj.colSize(idx) = 0;
            obj.colUp(idx) = 0;
            obj.colDown(idx) = 0;
            obj.rowRight(idx) = 0;
            obj.rowLeft(idx) = 0;
            obj.vplane(idx) = 0;
            obj.colCheck(idx) = 0;
            obj.XSC(idx) = 0;
            obj.YSC(idx) = 0;
            
            
        end
        %% Build new scattered data to pull missing objects from
        % This will use more sensitive object detection. *careful with false
        % postives!
        function [obj] = ReduceObjectDetectionThreshold(obj,image,raw)
            % First, identify the appropriate volume to search in image stack
            zmax = ceil(max(obj.Z(obj.plane==(1)))/raw.dataKey(10,1))+4;
            if zmax>size(image.MaskStack,3)
                zmax = size(image.RawStack,3);
            end
            zmin = floor(min(obj.Z(obj.plane==(1)))/raw.dataKey(10,1));
            d = 1;
            dz = 1;
            rp =  feature3dMB(image.MaskStack(:,:,zmin:zmax), [d d dz], [d d dz], [obj.s(1) obj.s(2) obj.s(3)],[1 0 0],2,0,0,0.1); %
            rp(:,1:2) = rp(:,1:2)*raw.dataKey(9,1);
            
            % remove detections near edges
            borderwidth = 2;% should be slightly larger than half the XY radius of an object
            tempr = rp;
            tempr = floor(tempr);
            removeidx = tempr(:,1)<=borderwidth | tempr(:,2)<=borderwidth | obj.s(1,2)-tempr(:,1)<=borderwidth | obj.s(2,2)-tempr(:,2)<=borderwidth;
            rp(removeidx,:) = [];
            
            %adjust X and Y data to match image orientation when
            %plotted (invert X, then switch X and Y)
            rp(:,1) = obj.s(1)*raw.dataKey(9,1)-rp(:,1);
            temp = rp(:,1);
            rp(:,1) = rp(:,2);
            rp(:,2) = temp;
            rp(:,3) = rp(:,3)*raw.dataKey(10,1)+(zmin-1)*raw.dataKey(10,1); %shift data back up
            
            % Perform a closest points search to eliminate duplicates
            [~,rpDist] = dsearchn(obj.r(obj.planeGroup==1,1:2),rp(:,1:2)); %xy check
            rp(rpDist<0.1,:) = [];
            
            [~,rpDist] = dsearchn(obj.r(obj.planeGroup==1,3),rp(:,3)); %z check
            rp(rpDist>0.5,:) = [];
            
            % Remove points which will not have potential matches within 3*
            % spacing
            
            tidx = obj.planeGroup ~= 1 & obj.colUp == 0;
            tr = obj.r(tidx,:);
            [~,rpDist] = dsearchn(tr(:,1:2),rp(:,1:2)); %suitale candidate check
            rp(rpDist>(raw.pSpaceXY*3),:) = [];
                        
            rpIdx = dsearchn(obj.r(:,1:3),rp(:,1:3));
            rpPlane = obj.plane(rpIdx);
            
            % Add new points to current lists
            rpOnes = ones(size(rp,1),1);
            rpZero = zeros(size(rp,1),1);
            rpZero3 = zeros(size(rp,1),3);
            rp(:,7:8) = 0;
            obj.r = cat(1,obj.r,rp);
            obj.X = cat(1,obj.X,rp(:,1));
            obj.Y = cat(1,obj.Y,rp(:,2));
            obj.Z = cat(1,obj.Z,rp(:,3));
            obj.rX = cat(1,obj.rX,rpZero);
            obj.rY = cat(1,obj.rY,rpZero);
            obj.rZ = cat(1,obj.rZ,rpZero);
            obj.dX = cat(1,obj.dX,rpZero);
            obj.dY = cat(1,obj.dY,rpZero);
            obj.dZ = cat(1,obj.dZ,rpZero);
            obj.dS = cat(1,obj.dS,rpZero);
            obj.Mass = cat(1,obj.Mass,rp(:,4));
            obj.MaxI = cat(1,obj.MaxI,rp(:,6));
            obj.rg = cat(1,obj.rg,rp(:,5));
            obj.PctAbove = cat(1,obj.PctAbove,rpZero);
            obj.row = cat(1,obj.row,rpZero);
            obj.col = cat(1,obj.col,rpZero);
            obj.colSize = cat(1,obj.colSize,rpZero);
            obj.colUp = cat(1,obj.colUp,rpZero);
            obj.colDown = cat(1,obj.colDown,rpZero);
            obj.rowRight = cat(1,obj.rowRight,rpZero);
            obj.rowLeft = cat(1,obj.rowLeft,rpZero);
            obj.plane = cat(1,obj.plane,rpPlane);
            obj.planeGroup = cat(1,obj.planeGroup,rpOnes);
            obj.vplane = cat(1,obj.vplane,rpZero);
            obj.colCheck = cat(1,obj.colCheck,rpZero);
            obj.Dl = cat(1,obj.Dl,rpZero);
            obj.NDl = cat(1,obj.NDl,rpZero);
            obj.rVa = cat(1,obj.rVa,rpZero3);
            obj.rVb = cat(1,obj.rVb,rpZero3);
            
        end
        %%
        function VerifyRows(obj,rows)
            rf = figure;
            hold on
            for i = 1:max(obj.row)
                rowIdx = obj.row==i;
                tempXY = zeros(nnz(rowIdx),3);
                tempXY(:,1) = obj.X(rowIdx);
                tempXY(:,2) = obj.Y(rowIdx);
                tempXY(:,3) = obj.Z(rowIdx);
                if abs(rows.V(1,1))<abs(rows.V(2,1))
                    [tempXY,~] = sortrows(tempXY,[2,1]);
                else
                    [tempXY,~] = sortrows(tempXY,[1,2]);
                end
                plot3(tempXY(:,1),tempXY(:,2),tempXY(:,3))
            end
            
            scatter3(obj.X(obj.ND),obj.Y(obj.ND),obj.Z(obj.ND),5,'white')
            scatter3(obj.X(obj.Dl),obj.Y(obj.Dl),obj.Z(obj.Dl),5,'black')
            scatter3(obj.X(obj.Dl&(obj.planeGroup==1)),obj.Y(obj.Dl&(obj.planeGroup==1)),obj.Z(obj.Dl&(obj.planeGroup==1)),5,'red')
            xlim([min(obj.X) ,max(obj.X)])
            ylim([min(obj.Y) ,max(obj.Y)])
            
            bcolor = 'black';
            fcolor = 'white';
            AxisFontSize = 12;
            LegendFontSize = 14;
            xt = 'Y \mum';% input('enter the xaxis label','s');
            yt = 'X \mum'; %input('enter the yaxis label','s');
            zt = 'Z \mum';
            label{1} = xlabel(xt);
            label{2} = ylabel(yt);
            label{3} = zlabel(zt);
            set(gca,'YMinorTick','on','color',bcolor)
            ytickformat('%.1f')
            le{1} = 'plane 1'; %input('enter the legend','s');
            le{2} = 'plane 2';
            ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,1,[0 0])
            %errorbar(meanDisplacements(1,1:3),meanDisplacements(2,1:3),'.','color',[0 0 0],'MarkerSize',1)
            
            legend off
            hold off
            view(-25,35)
            savefile = 'Rows.tif';
            export_fig(rf,savefile,'-native');
        end
        %%
        function viewColumn(obj)
            figure
            hold on
            for i = 1:max(obj.col)
                idx = obj.col == i;
                tempXY = zeros(nnz(idx),3);
                tempXY(:,1) = obj.X(idx);
                tempXY(:,2) = obj.Y(idx);
                tempXY(:,3) = obj.Z(idx);
                [tempXY,~] = sortrows(tempXY,3);
                plot3(tempXY(:,1),tempXY(:,2),tempXY(:,3))
            end
        end
    end
end
