classdef DispData3D
    properties
        ref
        refSC
        disp
        dispFilt
        dispMean
        MeanPlanes
        MeanTotal
        NegMax
        NegTotal
        PosMax
        PosTotal
        StdPlanes
        StdTotal
        XnoiseMean
        XnoiseStd
        XnoiseCO
        YnoiseMean
        YnoiseStd
        YnoiseCO
        SnoiseMean
        SnoiseStd
        SnoiseCO
        ZnoiseMean
        ZnoiseStd
        ZnoiseCO
        dispPF
        rowFits
        rowEnds
    end
    methods
        function obj = DispData3D
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Methods Functions 1,2 and 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %% Method 1 Fitting Non-Deformed Markers
        
        function obj = method1fit(obj,r,rows)
            rowFits = zeros(size(rows.m,1),3,2);
            progressbar('Method 1 Fits')
            
            %Generate linear fits of non-deformed markers in a single row
            for i = 1:size(rows.m,1)
                progressbar(i/(2*size(rows.m,1)))
                if nnz(intersect(rows.m(i,:),r.ND))>1
                    [A,B] = fitLine3D(r.X(intersect(rows.m(i,:),r.ND)),r.Y(intersect(rows.m(i,:),r.ND)),r.Z(intersect(rows.m(i,:),r.ND)));
                    rowFits(i,1:3,1) = A';
                    rowFits(i,1:3,2) = B';
                    
                elseif size(intersect(rows.m(i,:),r.ND),1)>0
                    [xyzFinal,rowP] = transLine3D(rows.V,r.r(intersect(rows.m(i,:),r.ND),1:3),r.s(1:3,1));
                    rowFits(i,1:3,1) = xyzFinal'-rows.V';
                    rowFits(i,1:3,2) = xyzFinal'+rows.V';
                else
                    [xyzFinal,rowP] = transLine3D(rows.V,r.r(rows.m(i,1:nnz(rows.m(i,:))),1:3),r.s(1:3,1));
                    rowFits(i,1:3,1) = xyzFinal'-rows.V';
                    rowFits(i,1:3,2) = xyzFinal'+rows.V';
                end
                
            end
            clear obj.rowFits
            
            %Extend the line segments generated above to beyond the image
            %frame. This is critical for properly finding the distance from
            %each line segment later.
            ty = (0 - rowFits(1,2,1))/(rowFits(1,2,2)-rowFits(1,2,1));
            tx = (0 - rowFits(1,1,1))/(rowFits(1,1,2)-rowFits(1,1,1));
            for i = 1:size(rowFits,1)
                progressbar(i+size(rows.m,1)/(2*size(rows.m,1)))
                if abs(ty)<abs(tx)
                    t = (0 - rowFits(i,2,1))/(rowFits(i,2,2)-rowFits(i,2,1));
                    obj.rowFits(i,1:3,1) = rowFits(i,1:3,1) + t*(rowFits(i,1:3,2)-rowFits(i,1:3,1));
                    t = (r.s(1,2) - rowFits(i,2,2))/(rowFits(i,2,1)-rowFits(i,2,2));
                    obj.rowFits(i,1:3,2) = rowFits(i,1:3,2) + t*(rowFits(i,1:3,1)-rowFits(i,1:3,2));
                else
                    t = (0 - rowFits(i,1,1))/(rowFits(i,1,2)-rowFits(i,1,1));
                    obj.rowFits(i,1:3,1) = rowFits(i,1:3,1) + t*(rowFits(i,1:3,2)-rowFits(i,1:3,1));
                    t = (r.s(2,2) - rowFits(i,1,2))/(rowFits(i,1,1)-rowFits(i,1,2));
                    obj.rowFits(i,1:3,2) = rowFits(i,1:3,2) + t*(rowFits(i,1:3,1)-rowFits(i,1:3,2));
                end
            end
            obj.rowFits = single(obj.rowFits);
        end
        
        
        %% Method 2 Translating Average Line
        function obj = method2fit(obj,m1,r,raw,im,rows)
            
            %% Revisiting Method One reinforced by Average slopes
            %Strategy: find closest and furthest points in row from one
            %point in rowFit and use those two points to determine if a
            %different fit using method 1 is necessary by querying whether
            %either end point falls within the list of non-deformed dots
            %under the cell.
            
            % Determine what the ends of each row are
            clear rowsEdgeDist obj.rowEnds
            progressbar('Method 2 Fits')
            rowS = size(rows.m,1);
            rowI = size(rows.Idx,1);
            ts = rowS + rowI + rowI + rowI + rowI + rowS;
            for i = 1:rowS
                progressbar(i/ts)
                clear distances currentRow currentMembers
                currentRow = rows.m(i,1:nnz(rows.m(i,:)));
                currentMembers = r.r(currentRow,1:3);
                for j = 1:size(currentMembers,1)
                    rowsEdgeDist(i,j) = pdist([currentMembers(j,1:2);m1.rowFits(i,1:2,1)]);
                end
                obj.rowEnds(i,1) = rows.m(i,find(rowsEdgeDist(i,:) == min(rowsEdgeDist(i,1:nnz(rowsEdgeDist(i,:))))));
                obj.rowEnds(i,2) = rows.m(i,find(rowsEdgeDist(i,:) == max(rowsEdgeDist(i,1:nnz(rowsEdgeDist(i,:))))));
                obj.rowEnds(i,3) = pdist([r.r(obj.rowEnds(i,1),1:3);r.r(obj.rowEnds(i,2),1:3)]);
            end
            
            % Determine which rows are useful in calculating an average row
            % slope Both ends must be outside of image.ADil's Cell and the
            % line width must be at least 70 percent of maximum line width
            for i = 1:rowI
                progressbar((i+rowS)/ts)
                rowPlanesMaxWidth(i,1) = max(obj.rowEnds(rows.Idx(i,1):rows.Idx(i,2),3));
                obj.rowEnds(rows.Idx(i,1):rows.Idx(i,2),5) = rowPlanesMaxWidth(i,1);
            end
            for i = 1:size(obj.rowEnds)
                progressbar((i+rowS+rowI)/ts)
                if r.NDl(obj.rowEnds(i,1))==1 && r.NDl(obj.rowEnds(i,2))==1 && obj.rowEnds(i,3)>obj.rowEnds(i,5)*.7 %|| obj.rowEnds(i,?)>20
                    obj.rowEnds(i,4) = 1;
                else
                    obj.rowEnds(i,4) = 0;
                end
            end
            
            % Calculate an average row slope per plane
            for i = 1:rowI
                progressbar((i+rowS+rowI+rowI)/ts)
                rowPlanesFits(i,1:3) = mean(m1.rowFits(find(obj.rowEnds(rows.Idx(i,1):rows.Idx(i,2),4)==1),1:3,2)-m1.rowFits(find(obj.rowEnds(rows.Idx(i,1):rows.Idx(i,2),4)==1),1:3,1));
                if size(find(obj.rowEnds(rows.Idx(i,1):rows.Idx(i,2),4)==1),1) <1
                    rowPlanesFits(i,1:3) = mean(m1.rowFits(find(obj.rowEnds(rows.Idx(i,1):rows.Idx(i,2),4)==0),1:3,2)-m1.rowFits(find(obj.rowEnds(rows.Idx(i,1):rows.Idx(i,2),4)==0),1:3,1));
                end
            end
            
            
            % View ends and average fits
            %             figure hold on
            %             scatter3(r.X(obj.rowEnds(:,1)),r.Y(obj.rowEnds(:,1)),r.Z(obj.rowEnds(:,1)))
            %             scatter3(r.X(obj.rowEnds(:,2)),r.Y(obj.rowEnds(:,2)),r.Z(obj.rowEnds(:,2)))
            %             scatter3(0,0,0) plot3([0 rowPlanesFits(1,1)],[0
            %             rowPlanesFits(1,2)],[0 rowPlanesFits(1,3)])
            %             axis([0 r.s(1,2) 0 r.s(2,2)])
            
            
            
            for i = 1:rowI
                progressbar((i+rowS+rowI+rowI+rowI)/ts)
                rowPlanesIdx2(i,1:pdist([rows.Idx(i,1);rows.Idx(i,2)])+1) = rows.Idx(i,1):1:rows.Idx(i,2);
            end
            %
            for i = 1:rowS
                progressbar((i+rowS+rowI+rowI+rowI+rowI)/ts)
                [rPIdxX rPIdxY] = find(rowPlanesIdx2 == i);
                if size(intersect(rows.m(i,:),r.ND),1)>0
                    %[xyzFinal,rowP] =
                    %transLine3D(rows.V,r.r(intersect(rows.m(i,:),r.ND),1:3),r.s(1:3,1));
                    [xyzFinal,rowP] = transLine3D(rowPlanesFits(rPIdxX,1:3)',r.r(intersect(rows.m(i,:),r.ND),1:3),r.s(1:3,1));
                    preRowFits(i,1:3,1) = xyzFinal-rowPlanesFits(rPIdxX,1:3)';
                    preRowFits(i,1:3,2) = xyzFinal+rowPlanesFits(rPIdxX,1:3)';
                else
                    [xyzFinal,rowP] = transLine3D(rowPlanesFits(rPIdxX,1:3)',r.r(obj.rowEnds(i,1:2),1:3),r.s(1:3,1));
                    preRowFits(i,1:3,1) = xyzFinal-rowPlanesFits(rPIdxX,1:3)';
                    preRowFits(i,1:3,2) = xyzFinal+rowPlanesFits(rPIdxX,1:3)';
                    
                end
                
            end
            % Extend the fit Lines passed the edge of the dot area - Method
            % 1 Can't believe this was a problem...
            
            for i = 1:size(preRowFits,1)
                if mean(abs(rowPlanesFits(:,1)))<mean(abs(rowPlanesFits(:,2)))
                    t = (0 - preRowFits(i,2,1))/(preRowFits(i,2,2)-preRowFits(i,2,1));
                    obj.rowFits(i,1:3,1) = preRowFits(i,1:3,1) + t*(preRowFits(i,1:3,2)-preRowFits(i,1:3,1));
                    t = (r.s(1,2) - preRowFits(i,2,2))/(preRowFits(i,2,1)-preRowFits(i,2,2));
                    obj.rowFits(i,1:3,2) = preRowFits(i,1:3,2) + t*(preRowFits(i,1:3,1)-preRowFits(i,1:3,2));
                else
                    t = (0 - preRowFits(i,1,1))/(preRowFits(i,1,2)-preRowFits(i,1,1));
                    obj.rowFits(i,1:3,1) = preRowFits(i,1:3,1) + t*(preRowFits(i,1:3,2)-preRowFits(i,1:3,1));
                    t = (r.s(2,2) - preRowFits(i,1,2))/(preRowFits(i,1,1)-preRowFits(i,1,2));
                    obj.rowFits(i,1:3,2) = preRowFits(i,1:3,2) + t*(preRowFits(i,1:3,1)-preRowFits(i,1:3,2));
                end
            end
            
        end
        
        
        %% Method 3 (1 and 2 combined)
        function obj = method3fit(obj,method1,method2,rows,raw)
            for i = 1:size(rows.m,1)
                if method2.rowEnds(i,4) == 1 && method2.rowEnds(i,3) > (raw.pSpaceXY * 10)
                    obj.rowFits(i,1:3,1:2) = method1.rowFits(i,1:3,1:2);
                else
                    obj.rowFits(i,1:3,1:2) = method2.rowFits(i,1:3,1:2);
                end
            end
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Displacement Functions 1 and 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         function obj = calcDisp(obj,r,rows)
%             obj.ref = zeros(r.l,3); for i = 1:size(obj.rowFits,1)
%                 if nnz(rows.NDm(i,:))>0
%                     n = nnz(rows.m(i,:)); temp =
%                     squeeze(obj.rowFits(i,:,:))'; [xy,distance,t_a] =
%                     distance2curve(temp,r.r(rows.m(i,1:n),1:3)); for j =
%                     1:size(xy,1)
%                         obj.ref(rows.m(i,j),1:3) = xy(j,1:3);
%                     end
%                 else
%                     for j = 1:nnz(rows.m(i,:))
%                         obj.ref(rows.m(i,j),1:3) = r.r(rows.m(i,j),1:3);
%                     end
%                 end
%             end
%         end
        
        function obj = calcDisp(obj,r,VPlanes)
            % XY component of reference Method 2
            for i = 1:r.l
                if r.row(i)>0
                    temp = squeeze(obj.rowFits(r.row(i),1:2,:))';
                    if r.rX(i)>0 %if a previous reference location exists:
                        %find the closest point on the row fit to the
                        %existing reference location
                        [xy,distance,t_a] = distance2curve(temp,[r.rX(i) r.rY(i)]);
                        obj.ref(i,1:2) = xy(1,1:2);
                    elseif r.colDown(i) >0
                        %otherwise, check to see if a reference location
                        %for the object below exists
                        if r.rX(r.colDown(i)) > 0
                            %if it does, use that reference location and
                            %the known
                            vppidx = max(r.planeGroup) - r.planeGroup(i);
                            tX = r.rX(r.colDown(i)) + VPlanes.vX(vppidx,1);
                            tY = r.rY(r.colDown(i)) + VPlanes.vY(vppidx,1);
                            [xy,distance,t_a] = distance2curve(temp,[tX tY]);
                            obj.ref(i,1:2) = xy(1,1:2);
                        else
                            obj.ref(i,1:2) = 0;
                        end
                    else
                        obj.ref(i,1:2) = 0;
                    end
                end
            end
            
            % Z component of reference Method 2
            for i = 1:r.l
                if r.row(i)>0
                differential = obj.rowFits(r.row(i),1:3,2)-obj.rowFits(r.row(i),1:3,1);
                differential2 = obj.rowFits(r.row(i),1:2,2) - obj.ref(i,1:2);
                differential3 = mean(differential(1,1:2)./differential2(1,1:2));
                obj.ref(i,3) = obj.rowFits(r.row(i),3,2)-differential(1,3)/differential3;
                end
            end
        end
        
        function [obj,r] = calcDispSC(obj,r,VPlanes)
            for i = 1:r.l
                if r.vplane(i)>0
                vppidx = max(r.planeGroup)-r.planeGroup(i);
                dx = sum(VPlanes.vX(r.vplane(i),1:vppidx));
                dy = sum(VPlanes.vY(r.vplane(i),1:vppidx));
                obj.refSC(i,1) = obj.ref(i,1)-dx;
                obj.refSC(i,2) = obj.ref(i,2)-dy;
                obj.refSC(i,3) = obj.ref(i,3);
                end
            end                        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Stats Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = dispStats(obj,plane,rows,r,planesLocFiltList)
            %% Calculate Average Displacement per Row per Plane Method 2
            for i = 1:size(rows.Idx,1)
                for j = 1:(rows.Idx(i,2)-rows.Idx(i,1))+1
                    obj.dispMean(j,i) = mean(abs(obj.disp(rows.NDm(j+(rows.Idx(i,1)-1),rows.NDm(j+(rows.Idx(i,1)-1),:)>0),3)));
                end
                obj.MeanPlanes(1,i) = mean(abs(obj.disp(intersect(r.ND,plane.final(1:nnz(plane.final(:,i)),i)),3)));
                obj.StdPlanes(1,i) = std(abs(obj.disp(intersect(r.ND,plane.final(1:nnz(plane.final(:,i)),i)),3)));
            end
            obj.MeanTotal= mean(abs(obj.disp(intersect(planesLocFiltList(:,1),r.ND),3)));
            obj.StdTotal = std(abs(obj.disp(intersect(planesLocFiltList(:,1),r.ND),3)));
        end
        
        function obj = dispNoise(obj,r,planesLocFiltList,image,scale)
            %% Filter out noise in Displacements Method 2
            %determine usable data
            useIdx = intersect(find(r.rX>0),intersect(planesLocFiltList(:,1),r.ND));
            
            obj.XnoiseMean = mean(abs(obj.disp(useIdx,1)));
            obj.XnoiseStd = std((obj.disp(useIdx,1)));
            obj.XnoiseCO = 2*obj.XnoiseStd;
            
            obj.YnoiseMean = mean(abs(obj.disp(useIdx,2)));
            obj.YnoiseStd = std((obj.disp(useIdx,2)));
            obj.YnoiseCO= 2*obj.YnoiseStd;
            
            obj.disp(:,4) = (obj.disp(:,1).^2 + obj.disp(:,2).^2).^(1/2);
            
            obj.SnoiseMean = mean(abs(obj.disp(useIdx,4)));
            obj.SnoiseStd = std((obj.disp(useIdx,4)));
            obj.SnoiseCO= 2*obj.SnoiseStd;
            
            obj.ZnoiseMean = mean(abs(obj.disp(useIdx,3)));
            obj.ZnoiseStd = std((obj.disp(useIdx,3)));
            obj.ZnoiseCO= 2*obj.ZnoiseStd;
            
            image.imgNBds = imageNoiseBounds(r.r,obj,image,scale);
            
            % Use noiseCutoff to filter data
            obj.dispFilt = obj.disp;
            for i = 1:size(r.r,1)
                try
                    
                    if abs(obj.dispFilt(i,3))<obj.ZnoiseCO && image.imgNBds(round(r.r(i,2)/scale),round(r.s(2,1)-r.r(i,1)/scale)) > 0
                        obj.dispFilt(i,3) = NaN;
                    elseif abs(obj.dispFilt(i,3))<obj.ZnoiseCO
                        obj.dispFilt(i,3) = 0;
                    end
                    
                catch
                    obj.dispFilt(i,3) = 0; %may need to correct this later, it will make some image border values zero
                end
            end
            obj.dispPF = find(abs(obj.disp(:,3))<obj.ZnoiseCO+obj.ZnoiseStd & abs(obj.disp(:,3))>obj.ZnoiseCO);
            clear rNbor
            radXY = 2.5;
            radZ = .75;
            for i = 1:size(r.r,1)
                topX = r.r(i,1)+ radXY;
                botX = r.r(i,1)- radXY;
                topY = r.r(i,2)+ radXY;
                botY = r.r(i,2)- radXY;
                topZ = r.r(i,3)+ radZ;
                botZ = r.r(i,3)- radZ;
                rNbor(i,1:size(find(r.r(:,1)<topX & r.r(:,1)>botX & r.r(:,2)<topY & r.r(:,2)>botY& r.r(:,3)<topZ & r.r(:,3)>botZ))) = find(r.r(:,1)<topX & r.r(:,1)>botX & r.r(:,2)<topY & r.r(:,2)>botY& r.r(:,3)<topZ & r.r(:,3)>botZ);
            end
            
            
            for i = 1:size(obj.dispPF,1)
                A = obj.dispPF(i,1);
                obj.dispPF(i,2)= obj.disp(obj.dispPF(i,1),3);
                obj.dispPF(i,3) = mean(obj.dispFilt(rNbor(A,rNbor(A,:)~=0&(rNbor(A,:)~=A)),3));
                obj.dispPF(i,4) = std(obj.dispFilt(rNbor(A,rNbor(A,:)~=0&(rNbor(A,:)~=A)),3));
                %filt 2 std from mean
                if pdist([obj.dispPF(i,2);obj.dispPF(i,3)]) > ((obj.dispPF(i,4)))
                    obj.dispPF(i,5) = 0;
                else
                    obj.dispPF(i,5) = 1;
                end
                
                %     if obj.dispPF(i,5) == 0
                %         obj.dispFilt(obj.dispPF(i,1),3) = 0;
                %     end
%                 
                
            end
        end
        
        function printStats(obj,planesGroups,planesLocTxt,planesLoc2,method)
            
            planesLocComma = fopen(strcat('Planes Data ',method,'.txt'),'wt');
            p1Format = 'Plane\tDistance\tMean\tStd\tPositive Total\tPositive Max\tNegative Total\tNegative Max \n';
            p2Format = '%1.0f\t%.2f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f \n';
            fprintf(planesLocTxt,p1Format);
            for i = 1:size(planesGroups,1)
                fprintf(planesLocTxt,p2Format,i,mean(planesLoc2(1,planesGroups(i,(planesGroups(i,:)>0)))),...
                    mean(obj.MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))),mean(obj.StdPlanes(1,planesGroups(i,(planesGroups(i,:)>0))))...
                    ,obj.PosTotal(i,1),obj.PosMax(i,1),obj.NegTotal(i,1),obj.NegMax(i,1));
            end
            fclose(planesLocComma);
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Plotting Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ViewMethodRef(obj,r,method)
            AxisFontSize = 24;
            AxisTitleFontSize = 24;
            LegendFontSize = 14;
            colOptions{1,1} = 'white';
            colOptions{2,1} = 'black';
            colOptions{1,2} = 'black';
            colOptions{2,2} = 'white';
            
            for k = 1:size(colOptions,2)
                fcolor = colOptions{1,k};
                bcolor = colOptions{2,k};
                refLocs = figure;
                set(gcf,'unit','pixels','position',[500,200,1000,500])
                
                hold on
                scatter3(obj.ref(:,1),obj.ref(:,2),obj.ref(:,3))
                for i = 1:size(obj.rowFits,1)
                    plot3([obj.rowFits(i,1,1) obj.rowFits(i,1,2)],[obj.rowFits(i,2,1) obj.rowFits(i,2,2)],[obj.rowFits(i,3,1) obj.rowFits(i,3,2)])
                end
                axis([0 r.s(1,2) 0 r.s(2,2)])
                xt = 'X-axis (\mum)';% input('enter the xaxis label','s');
                yt = 'Y-axis (\mum)'; %input('enter the yaxis label','s');
                zt = 'Z-axis (\mum)'; %input('enter the yaxis label','s');
                le{1} = 'Ref. Locs'; %input('enter the legend','s');
                %             le{2} = 'Time-Lapse 1'; le{3} = 'Time-Lapse
                %             2';
                label{1} = xlabel(xt);
                label{2} = ylabel(yt);
                label{3} = zlabel(zt);
                view(1) = -45; view(2) = 8;
                ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,1,view)
                
                %Export Image
                mkdir '3D Plots'
                title = ['\' method ' Reference Locations ' fcolor ' on ' bcolor];
                savefile = [cd '\3D Plots' title];
                export_fig(refLocs,savefile,'-native');
            end
            
        end
        %%
        function ViewRowFits(obj,r,rows,method)
            AxisFontSize = 24;
            AxisTitleFontSize = 24;
            LegendFontSize = 14;
            colOptions{1,1} = 'white';
            colOptions{2,1} = 'black';
            colOptions{1,2} = 'black';
            colOptions{2,2} = 'white';
            
            for k = 1:size(colOptions,2)
                fcolor = colOptions{1,k};
                bcolor = colOptions{2,k};
                Locs = figure;
                set(gcf,'unit','pixels','position',[500,200,1000,500])                               
                
                hold on
                
                for i = 1:max(r.plane)
                    idx = r.plane == i;
                    scatter3(r.X(idx),r.Y(idx),r.Z(idx),2)
                end
                
                for i = 1:size(obj.rowFits,1)
                    plot3([obj.rowFits(i,1,1),obj.rowFits(i,1,2)],[obj.rowFits(i,2,1),obj.rowFits(i,2,2)],[obj.rowFits(i,3,1),obj.rowFits(i,3,2)])
                end
                
                
                axis([0 r.s(2,2) 0 r.s(1,2)])
                zlim([min(r.Z)-1 max(r.Z)+1])
                xlim([19 35])
                ylim([0 80])
                
                xt = 'Y (\mum)';% input('enter the xaxis label','s');
                yt = 'X (\mum)'; %input('enter the yaxis label','s');
                zt = 'Z (\mum)'; %input('enter the yaxis label','s');
                le{1} = 'Centroids'; %input('enter the legend','s');
                %             le{2} = 'Time-Lapse 1'; le{3} = 'Time-Lapse
                %             2';
                label{1} = xlabel(xt);
                label{2} = ylabel(yt);
                label{3} = zlabel(zt);
                view(1) = -9; view(2) = 24;
                ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,1,view)
                
                %Export Image
                mkdir '3D Plots'
                title = ['\' method ' Marker Locations ' fcolor ' on ' bcolor];
                savefile = [cd '\3D Plots' title];
                export_fig(Locs,savefile,'-native');
            end
        end
        %%
        function ViewQuiverPlot(obj,r,method)
            AxisFontSize = 24;
            AxisTitleFontSize = 24;
            LegendFontSize = 14;
            colOptions{1,1} = 'white';
            colOptions{2,1} = 'black';
            colOptions{1,2} = 'black';
            colOptions{2,2} = 'white';
            
            for k = 1:size(colOptions,2)
                fcolor = colOptions{1,k};
                bcolor = colOptions{2,k};
                Quiver = figure;
                set(gcf,'unit','pixels','position',[500,200,1000,500])
                hold on
                idx = obj.ref(:,1) ~= 0;
                quiver3(obj.ref(idx,1),obj.ref(idx,2),obj.ref(idx,3),obj.disp(idx,1),obj.disp(idx,2),obj.disp(idx,3),0)
                axis([0 r.s(2,2) 0 r.s(1,2)])
                xt = 'Y (\mum)';% input('enter the xaxis label','s');
                yt = 'X (\mum)'; %input('enter the yaxis label','s');
                zt = 'Z (\mum)'; %input('enter the yaxis label','s');
                le{1} = 'Disp. Vector'; %input('enter the legend','s');
                %             le{2} = 'Time-Lapse 1'; le{3} = 'Time-Lapse
                %             2';
                label{1} = xlabel(xt);
                label{2} = ylabel(yt);
                label{3} = zlabel(zt);
                
                    zlim([min(r.Z)-1 max(r.Z)+1])
                xlim([19 35])
                ylim([0 80])
                view(1) = -9; view(2) = 24;
                ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,1,view)
                
                %Export Image
                mkdir '3D Plots'
                title = ['\' method ' Quiver Plot ' fcolor ' on ' bcolor];
                savefile = [cd '\3D Plots' title];
                export_fig(Quiver,savefile,'-native');
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Histogram Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function NoiseHists(obj,planesLocFiltList,r,method)
            mkdir('Histograms')
            %% Create Histogram
            AxisFontSize = 28;
            AxisTitleFontSize = 28;
            LegendFontSize = 20;
            colOptions{1,1} = 'white';
            colOptions{2,1} = 'black';
            colOptions{1,2} = 'black';
            colOptions{2,2} = 'white';
            useIdx = intersect(find(r.rX>0),intersect(planesLocFiltList(:,1),r.ND));
            
            for k = 1:size(colOptions,2)
                fcolor = colOptions{1,k};
                bcolor = colOptions{2,k};
                zdispdist = figure;
                hold on
                bins = [-400:20:400];
                SigColor = [.1 .7 .7];
                histogram(obj.disp(useIdx,3)*1000,bins,'FaceColor',fcolor,'Normalization','probability','EdgeColor',fcolor)
                histmax = max(histcounts(obj.disp(useIdx,3)*1000,bins,'Normalization','probability'));
                p1 = plot([obj.ZnoiseCO/2*1000 obj.ZnoiseCO/2*1000],[0 .5],'color',[.7 .7 .7],'linestyle','--','linewidth',1);
                plot([(obj.ZnoiseCO*-1)/2*1000 (obj.ZnoiseCO*-1)/2*1000],[0 .5],'color',[.7 .7 .7],'linestyle','--','linewidth',1)
                p2 = plot([obj.ZnoiseCO*1000 obj.ZnoiseCO*1000],[0 .5],'color',fcolor,'linestyle','--','linewidth',1);
                plot([obj.ZnoiseCO*-1*1000 obj.ZnoiseCO*-1*1000],[0 .5],'color',fcolor,'linestyle','--','linewidth',1)
                set(gca,'fontsize',28,'YMinorTick','on')
                xt = 'Z-Displacement (nm)';% input('enter the xaxis label','s');
                yt = 'Probability'; %input('enter the yaxis label','s');
                tt = 'Line-Profile Displacements';%input('enter the title','s');
                le{1} = {'0'};
                le2 = ['\color{' fcolor '}\sigma']; %input('enter the legend','s');
                le3 = ['\color{' fcolor '}2*\sigma'];
                label{1} = xlabel(xt);
                label{2} = ylabel(yt);
                %tl = title(tt);
                ytickformat('%.2f')
                ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,0,360)
                leg = legend([p1 p2],le2,le3,'location',[ .8 .8 .1 .1]);
                leg.EdgeColor = fcolor;
                leg.FontSize = LegendFontSize;
                axis([-400 400 0 .3])
                
                text((double(obj.ZnoiseCO)*1000+20),(.15),strcat('2\sigma= ',num2str(round(obj.ZnoiseCO*1000,0)),'nm'),'color',fcolor ,'fontsize',18)
                text((double(obj.ZnoiseCO)*1000+20),(.12),'Noise Cutoff','color',fcolor ,'fontsize',16 )
                legend boxoff
                title = ['Z-Displacement Histogram '  method ' ' fcolor ' on ' bcolor];
                filePath = cd;
                savefile = [filePath '\Histograms\' title];
                export_fig(zdispdist,savefile,'-native');
            end
        end
    end
end
