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
        noiseMean
        noiseStd
        noiseCutoff
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
        
        function obj = method1fit(obj,r,rows,rowV)
            rowFits = zeros(size(rows,1),3,2);
            for i = 1:size(rows,1)
                if size(intersect(rows(i,:),r.ND),1)>1
                    [A,B] = fitLine3D(r.X(intersect(rows(i,:),r.ND)),r.Y(intersect(rows(i,:),r.ND)),r.Z(intersect(rows(i,:),r.ND)));
                    rowFits(i,1:3,1) = A';
                    rowFits(i,1:3,2) = B';
                    
                elseif size(intersect(rows(i,:),r.ND),1)>0
                    [xyzFinal,rowP] = transLine3D(rowV,r.r(intersect(rows(i,:),r.ND),1:3),r.s(1:3,1));
                    rowFits(i,1:3,1) = xyzFinal'-rowV;
                    rowFits(i,1:3,2) = xyzFinal'+rowV;
                else
                    [xyzFinal,rowP] = transLine3D(rowV,r.r(rows(i,1:nnz(rows(i,:))),1:3),r.s(1:3,1));
                    rowFits(i,1:3,1) = xyzFinal'-rowV;
                    rowFits(i,1:3,2) = xyzFinal'+rowV;
                end
            end
            clear obj.rowFits
            
            ty = (0 - rowFits(1,2,1))/(rowFits(1,2,2)-rowFits(1,2,1));
            tx = (0 - rowFits(1,1,1))/(rowFits(1,1,2)-rowFits(1,1,1));
            for i = 1:size(rowFits,1)
                if abs(ty)<abs(tx)
                    t = (0 - rowFits(i,2,1))/(rowFits(i,2,2)-rowFits(i,2,1));
                    obj.rowFits(i,1:3,1) = rowFits(i,1:3,1) + t*(rowFits(i,1:3,2)-rowFits(i,1:3,1));
                    t = (r.s(2,2) - rowFits(i,2,2))/(rowFits(i,2,1)-rowFits(i,2,2));
                    obj.rowFits(i,1:3,2) = rowFits(i,1:3,2) + t*(rowFits(i,1:3,1)-rowFits(i,1:3,2));
                else
                    t = (0 - rowFits(i,1,1))/(rowFits(i,1,2)-rowFits(i,1,1));
                    obj.rowFits(i,1:3,1) = rowFits(i,1:3,1) + t*(rowFits(i,1:3,2)-rowFits(i,1:3,1));
                    t = (r.s(2,1) - rowFits(i,1,2))/(rowFits(i,1,1)-rowFits(i,1,2));
                    obj.rowFits(i,1:3,2) = rowFits(i,1:3,2) + t*(rowFits(i,1:3,1)-rowFits(i,1:3,2));
                end
            end
            obj.rowFits = single(obj.rowFits);
        end
        
        
        %% Method 2 Translating Average Line
        function obj = method2fit(obj,m1,r,raw,im,rows,rowPlanesIdx)
            
            %% Revisiting Method One reinforced by Average slopes
            %Strategy: find closest and furthest points in row from one point in rowFit
            %and use those two points to determine if a different fit using method 1 is
            %necessary by querying whether either end point falls within the list of
            %non-deformed dots under the cell.
            
            % Determine what the ends of each row are
            clear rowsEdgeDist obj.rowEnds
            for i = 1:size(rows,1)
                clear distances currentRow currentMembers
                currentRow = rows(i,1:nnz(rows(i,:)));
                currentMembers = r.r(currentRow,1:3);
                for j = 1:size(currentMembers,1)
                    rowsEdgeDist(i,j) = pdist([currentMembers(j,1:2);m1.rowFits(i,1:2,1)]);
                end
                obj.rowEnds(i,1) = rows(i,find(rowsEdgeDist(i,:) == min(rowsEdgeDist(i,1:nnz(rowsEdgeDist(i,:))))));
                obj.rowEnds(i,2) = rows(i,find(rowsEdgeDist(i,:) == max(rowsEdgeDist(i,1:nnz(rowsEdgeDist(i,:))))));
                obj.rowEnds(i,3) = pdist([r.r(obj.rowEnds(i,1),1:3);r.r(obj.rowEnds(i,2),1:3)]);
            end
            
            % Determine which rows are useful in calculating an average row slope
            % Both ends must be outside of image.ADil's Cell and the line width
            % must be at least 70 percent of maximum line width
            for i = 1:size(rowPlanesIdx,1)
                rowPlanesMaxWidth(i,1) = max(obj.rowEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),3));
                obj.rowEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),5) = rowPlanesMaxWidth(i,1);
            end
            for i = 1:size(obj.rowEnds)
                
                if im(ceil(r.Y(obj.rowEnds(i,1))/raw.dataKey(9,1)),ceil(r.X(obj.rowEnds(i,1))/raw.dataKey(9,1)))>0 && im(ceil(r.Y(obj.rowEnds(i,2))/raw.dataKey(9,1)),ceil(r.X(obj.rowEnds(i,2))/raw.dataKey(9,1)))>0 && obj.rowEnds(i,3)>obj.rowEnds(i,5)*.7 %|| obj.rowEnds(i,?)>20
                    obj.rowEnds(i,4) = 1;
                else
                    obj.rowEnds(i,4) = 0;
                end
            end
            
            % Calculate an average row slope per plane
            for i = 1:size(rowPlanesIdx,1)
                rowPlanesFits(i,1:3) = mean(m1.rowFits(find(obj.rowEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),4)==1),1:3,2)-m1.rowFits(find(obj.rowEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),4)==1),1:3,1));
                if size(find(obj.rowEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),4)==1),1) <1
                    rowPlanesFits(i,1:3) = mean(m1.rowFits(find(obj.rowEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),4)==0),1:3,2)-m1.rowFits(find(obj.rowEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),4)==0),1:3,1));
                end
            end
            
            
            % View ends and average fits
%             figure
%             hold on
%             scatter3(r.X(obj.rowEnds(:,1)),r.Y(obj.rowEnds(:,1)),r.Z(obj.rowEnds(:,1)))
%             scatter3(r.X(obj.rowEnds(:,2)),r.Y(obj.rowEnds(:,2)),r.Z(obj.rowEnds(:,2)))
%             scatter3(0,0,0)
%             plot3([0 rowPlanesFits(1,1)],[0 rowPlanesFits(1,2)],[0 rowPlanesFits(1,3)])
%             axis([0 r.s(1,2) 0 r.s(2,2)])
            
            
            
            for i = 1:size(rowPlanesIdx,1)
                rowPlanesIdx2(i,1:pdist([rowPlanesIdx(i,1);rowPlanesIdx(i,2)])+1) = rowPlanesIdx(i,1):1:rowPlanesIdx(i,2);
            end
            %
            for i = 1:size(rows,1)
                [rPIdxX rPIdxY] = find(rowPlanesIdx2 == i);
                if size(intersect(rows(i,:),r.ND),1)>0
                    %[xyzFinal,rowP] = transLine3D(rowV,r(rowsNDCU(i,1:n),1:3),imageSize);
                    [xyzFinal,rowP] = transLine3D(rowPlanesFits(rPIdxX,1:3),r.r(intersect(rows(i,:),r.ND),1:3),r.s(1:3,1));
                    preRowFits(i,1:3,1) = xyzFinal-rowPlanesFits(rPIdxX,1:3)';
                    preRowFits(i,1:3,2) = xyzFinal+rowPlanesFits(rPIdxX,1:3)';
                else
                    [xyzFinal,rowP] = transLine3D(rowPlanesFits(rPIdxX,1:3),r.r(obj.rowEnds(i,1:2),1:3),r.s(1:3,1));
                    preRowFits(i,1:3,1) = xyzFinal-rowPlanesFits(rPIdxX,1:3)';
                    preRowFits(i,1:3,2) = xyzFinal+rowPlanesFits(rPIdxX,1:3)';
                    i
                end
                
            end
            % Extend the fit Lines passed the edge of the dot area - Method 1
            % Can't believe this was a problem...
            
            for i = 1:size(preRowFits,1)
                if mean(abs(rowPlanesFits(:,1)))<mean(abs(rowPlanesFits(:,2)))
                    t = (0 - preRowFits(i,2,1))/(preRowFits(i,2,2)-preRowFits(i,2,1));
                    obj.rowFits(i,1:3,1) = preRowFits(i,1:3,1) + t*(preRowFits(i,1:3,2)-preRowFits(i,1:3,1));
                    t = (r.s(2,2) - preRowFits(i,2,2))/(preRowFits(i,2,1)-preRowFits(i,2,2));
                    obj.rowFits(i,1:3,2) = preRowFits(i,1:3,2) + t*(preRowFits(i,1:3,1)-preRowFits(i,1:3,2));
                else
                    t = (0 - preRowFits(i,1,1))/(preRowFits(i,1,2)-preRowFits(i,1,1));
                    obj.rowFits(i,1:3,1) = preRowFits(i,1:3,1) + t*(preRowFits(i,1:3,2)-preRowFits(i,1:3,1));
                    t = (r.s(1,2) - preRowFits(i,1,2))/(preRowFits(i,1,1)-preRowFits(i,1,2));
                    obj.rowFits(i,1:3,2) = preRowFits(i,1:3,2) + t*(preRowFits(i,1:3,1)-preRowFits(i,1:3,2));
                end
            end
            
        end
        
        
        %% Method 3 (1 and 2 combined)
        function obj = method3fit(obj,method1,method2,rows)
            for i = 1:size(rows,1)
                if method2.rowEnds(i,4) == 1
                    obj.rowFits(i,1:3,1:2) = method1.rowFits(i,1:3,1:2);
                else
                    obj.rowFits(i,1:3,1:2) = method2.rowFits(i,1:3,1:2);
                end
            end
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Displacement Functions 1 and 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = calcDisp(obj,r,rows,rowsNDCU)
            obj.ref = zeros(r.l,3);
            for i = 1:size(obj.rowFits,1)
                if nnz(rowsNDCU(i,:))>0
                    n = nnz(rows(i,:));
                    temp = squeeze(obj.rowFits(i,:,:))';
                    [xy,distance,t_a] = distance2curve(temp,r.r(rows(i,1:n),1:3));
                    for j = 1:size(xy,1)
                        obj.ref(rows(i,j),1:3) = xy(j,1:3);
                    end
                else
                    for j = 1:nnz(rows(i,:))
                        obj.ref(rows(i,j),1:3) = r.r(rows(i,j),1:3);
                    end
                end
            end
        end
        
        function obj = calcDispSC(obj,r,shear)
            % XY component of reference Method 2
            for i = 1:r.l
                temp = squeeze(obj.rowFits(r.row(i),1:2,:))';
                if r.col(i)>0
                    [xy,distance,t_a] = distance2curve(temp,[shear.rawX1(r.col(i)) shear.rawY1(r.col(i))]);
                    obj.refSC(i,1:2) = xy(1,1:2) + [r.XSC(i) r.YSC(i)];
                else
                    obj.refSC(i,1:2) = obj.ref(i,1:2);
                end
            end
            
            % Z component of reference Method 2
            for i = 1:r.l
                differential = obj.rowFits(r.row(i),1:3,2)-obj.rowFits(r.row(i),1:3,1);
                differential2 = obj.rowFits(r.row(i),1:2,2) - obj.refSC(i,1:2);
                differential3 = mean(differential(1,1:2)./differential2(1,1:2));
                obj.refSC(i,3) = obj.rowFits(r.row(i),3,2)-differential(1,3)/differential3;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Stats Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = dispStats(obj,plane,rowPlanesIdx,r,rowsNDCU,planesLocFiltList)
            [obj.dispMean,obj.MeanPlanes,obj.StdPlanes,obj.MeanTotal,obj.StdTotal] = quickRowStats(obj.disp,plane.final,rowPlanesIdx,r.ND,rowsNDCU,planesLocFiltList);
        end
        
        function obj = dispNoise(obj,r,planesLocFiltList)
            [obj.noiseMean,obj.noiseStd,obj.noiseCutoff,obj.dispFilt,obj.dispPF] = rowNoiseCalc(r.r,obj.disp,planesLocFiltList,r.ND);
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
        function ViewMethodRef(obj,r)
            figure
            hold on
            scatter3(obj.refSC(:,1),obj.refSC(:,2),obj.refSC(:,3))
            for i = 1:size(obj.rowFits,1)
                plot3([obj.rowFits(i,1,1) obj.rowFits(i,1,2)],[obj.rowFits(i,2,1) obj.rowFits(i,2,2)],[obj.rowFits(i,3,1) obj.rowFits(i,3,2)])
            end
            axis([0 r.s(1,2) 0 r.s(2,2)])
        end
        
        function ViewRowFits(obj,r,rowPlanes)
            figure
            hold on
            for j = 1:size(rowPlanes,3)
                for i = 1:size(rowPlanes,1)
                    n = nnz(rowPlanes(i,:,j));
                    scatter3(r.X(rowPlanes(i,1:n,j)),r.Y(rowPlanes(i,1:n,j)),r.Z(rowPlanes(i,1:n,j)))
                end
            end
            
            for i = 1:size(obj.rowFits,1)
                plot3([obj.rowFits(i,1,1) obj.rowFits(i,1,2)],[obj.rowFits(i,2,1) obj.rowFits(i,2,2)],[obj.rowFits(i,3,1) obj.rowFits(i,3,2)])
            end
            scatter3(0,0,0)
            axis([0 r.s(1,2) 0 r.s(2,2)])
        end
        
        function ViewQuiverPlot(obj,r)
            figure
            quiver3(obj.refSC(:,1),r.s(2,2)-obj.refSC(:,2),obj.refSC(:,3),obj.disp(:,1),obj.disp(:,2)*-1,obj.disp(:,3),0)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Histogram Functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function NoiseHists(obj,planesLocFiltList,r,method)
            %% Create Histogram
            zdispdist = figure;
            hold on
            bins = [-400:20:400];
            SigColor = [.1 .7 .7];
            histogram(obj.disp(intersect(planesLocFiltList(:,1),r.ND),3)*1000,bins,'FaceColor',[.6 .6 .6],'Normalization','probability')
            histmax = max(histcounts(obj.disp(intersect(planesLocFiltList(:,1),r.ND),3)*1000,bins,'Normalization','probability'));
            p1 = plot([obj.noiseCutoff/2*1000 obj.noiseCutoff/2*1000],[0 round(histmax,2)+.01],'color',[.3 .3 .3],'linestyle','--','linewidth',1);
            plot([(obj.noiseCutoff*-1)/2*1000 (obj.noiseCutoff*-1)/2*1000],[0 round(histmax,2)+.01],'color',[.3 .3 .3],'linestyle','--','linewidth',1)
            p2 = plot([obj.noiseCutoff*1000 obj.noiseCutoff*1000],[0 round(histmax,2)+.01],'color',SigColor ,'linestyle','--','linewidth',1);
            plot([obj.noiseCutoff*-1*1000 obj.noiseCutoff*-1*1000],[0 round(histmax,2)+.01],'color',SigColor ,'linestyle','--','linewidth',1)
            set(gca,'fontsize',28)
            xt = 'Z-Displacement (nm)';% input('enter the xaxis label','s');
            yt = 'Probability'; %input('enter the yaxis label','s');
            tt = 'Line-Profile Displacements';%input('enter the title','s');
            le = '\sigma'; %input('enter the legend','s');
            le2 = '2*\sigma';
            le3 = 'Cell Border';
            xl = xlabel(xt);
            yl = ylabel(yt);
            %tl = title(tt);
            
            set(xl, 'fontweight','bold','fontsize',28);
            set(yl,'fontweight','bold','fontsize',28);
            leg = legend([p1 p2],le,le2,'location','northwest');
            leg.FontSize = 20;
            axis([-400 400 0 round(histmax,2)+.01])
            
            text((double(obj.noiseCutoff)*1.01*1000),(histmax*.5),strcat('2\sigma= ',num2str(round(obj.noiseCutoff*1000,0)),'nm'),'color',SigColor ,'fontsize',20)
            text((double(obj.noiseCutoff)*1.01*1000),(histmax*.5-.01),'Noise Cutoff','color',SigColor ,'fontsize',18 )
            
            title = strcat('Z-Displacement Histogram',method);
            filePath = cd;
            savefile = [filePath '\Histograms\' title];
            export_fig(zdispdist,savefile,'-native');
            
        end
    end
end
