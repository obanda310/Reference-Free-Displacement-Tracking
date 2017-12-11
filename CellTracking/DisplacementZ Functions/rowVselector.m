function [rowV, rROI, rROIplane] = rowVselector(r,raw,image)
%% Dots in Cropped Region
    clear rNDB neighbors
    %ND is non-deformed
    radXY =  2.5; %microns
    radZ = .3;
    rROIBounds = image.SquareBounds*raw.dataKey(9,1);
    rROIBounds(1,5:6) = (rROIBounds(1,1:2) + rROIBounds(1,3:4))/2;
    rROIImage = image.TDots(image.SquareBounds(1,2):image.SquareBounds(1,4),image.SquareBounds(1,1):image.SquareBounds(1,3));
    rROIInd(:,1)  = (r.X(:)>rROIBounds(1,1) & r.X(:)<rROIBounds(1,3));
    rROIInd(:,2)  = (r.Y(:)>rROIBounds(1,2) & r.Y(:)<rROIBounds(1,4));
    rROIInd(:,3) = rROIInd(:,1) .* rROIInd(:,2);
    rROIIndFinal = find(rROIInd(:,3));
    rND = (r.r(rROIIndFinal,1:end));
    rROI = RawData3D(rROIImage,raw,rND);
    rROI = TranscribeR(rROI);
    %viewDetections(rROI,raw,rROIBounds,rROIImage);
    %% Find plane with least fluctuations
    
    rROIplane = PlanesData(rROI,radXY,radZ);
    rROIplane = growPlanes(rROIplane,rROI);
    
    %%  
    %Determine best (most populated) plane for an approximate non-deformed region    
    for i = 1:size(rROIplane.raw,2)
    rROIpNNZ(i) = nnz(rROIplane.raw(:,i));
    rROIpB = find(rROIpNNZ==max(rROIpNNZ),1,'first'); %find most populated
    end
    rNDB = rROI.r(rROIplane.raw(:,rROIpB),:);
    
    %Find Center dot (likely to have 4 equidistant neighbors)
    
        differences(:,1) = rNDB(:,1)-rROIBounds(1,5);
        differences(:,2) = rNDB(:,2)-rROIBounds(1,6);
        differences(:,4) = sqrt(differences(:,1).^2 + differences(:,2).^2);
    
    best = find(differences(:,4)==min(differences(:,4)));
    %%
    % figure
%     scatter3(rNDB(:,1),rNDB(:,2),rNDB(:,3))
%     hold on
%     scatter3(rNDB(best,1),rNDB(best,2),rNDB(best,3))
%     hold off
    %%
    k = 0;
    count = 0;
    while k == 0
        %Find the 4 'equidistant' neighbors
        clear differences sortedNew sortedOrig neighbors
        for i = 1:size(rNDB,1)
            differences(i,1:3) = rNDB(i,1:3)-rNDB(best,1:3);
            differences(i,4) = sqrt(differences(i,1)^2 + differences(i,2)^2);
        end
        [sortedNew, sortedOrig] = sort(differences(:,4));
        if std(sortedNew(2:5))<.5
            neighbors = rNDB(sortedOrig(2:5),:);
            neighbors(1:4,9) = sortedOrig(2:5);
            k = 1;
        elseif count > 20
            k = 1;
            disp('Could not find a suitable candidate for line fit')
        else
            newGuess = 2+round(3*rand());
            best = sortedOrig(newGuess);
            count = count +1;
        end
        
    end
    %%
    % figure
    % scatter3(rNDB(:,1),rNDB(:,2),rNDB(:,3))
    % hold on
    % scatter3(neighbors(:,1),neighbors(:,2),neighbors(:,3))
    % scatter3(rNDB(best,1),rNDB(best,2),rNDB(best,3))
    % hold off
    % % trckText = strcat('\leftarrow ',trckNum);
    % % text(lub(nghbrs(i,2),1)-(cntrPt(1,1)-fSizeXmin),lub(nghbrs(i,2),2)-(cntrPt(1,2)-fSizeYmin),lub(nghbrs(i,2),6),trckText,'Color','red')
    %
    %%
%     m = 15;
%     figure
%     scatter3(rNDB(:,1),rNDB(:,2),rNDB(:,3))
%     hold on
%     scatter3(rNDB(sortedOrig(1:m,1),1),rNDB(sortedOrig(1:m,1),2),rNDB(sortedOrig(1:m,1),3))
%     m=10;
%     scatter3(rNDB(sortedOrig(1:m,1),1),rNDB(sortedOrig(1:m,1),2),rNDB(sortedOrig(1:m,1),3))
%     m=5;
%     scatter3(rNDB(sortedOrig(1:m,1),1),rNDB(sortedOrig(1:m,1),2),rNDB(sortedOrig(1:m,1),3))
%     scatter(0,0,0)
%      hold off
    %%
    % Pair off neighbors
    clear differences
    for i = 1:4
        for j = 1:4
            differences(j,1:3) = neighbors(i,1:3)-neighbors(j,1:3);
            differences(j,4) = sqrt(differences(j,1)^2 + differences(j,2)^2 + differences(j,3)^2);
        end
        neighbors(i,10) = find(differences(:,4)==max(differences(:,4)));
    end
    
    used = zeros(1,4);
    for i = 1:4
        if ismember(i,used)
        else
            dFit{i}(1,1:3) = neighbors(i,1:3);
            dFit{i}(2,1:3) = neighbors(neighbors(i,10),1:3);
            dFit{i}(3,1:3) = rNDB(best,1:3);
            used(i,1) = neighbors(i,10);
        end
    end
    
    %%
    
    v1 = (dFit{1}(1,1:2) - dFit{1}(2,1:2))/2;
    v2 = (dFit{2}(1,1:2) - dFit{2}(2,1:2))/2;
    %%
    clear v1row
    v1row = best;
    for i = 1:2
        count = 1;
        dv1 = rNDB(best,1:2);
        while (dv1(1,1) < rROIBounds(1,3) && dv1(1,1) > rROIBounds(1,1) && dv1(1,2) < rROIBounds(1,4) && dv1(1,2) > rROIBounds(1,2)) == 1       
                clear differences
                dv1 = dv1 + v1*(-1^i);
                dv11(1:size(rNDB,1),1)=dv1(1,1);
                dv11(1:size(rNDB,1),2)=dv1(1,2);
                differences(:,1:2) = (rNDB(:,1:2) - dv11(:,1:2));
                differences(:,3) = sqrt(differences(:,1).^2+differences(:,2).^2);
                [dvSortNew,dvSortOrig] = sort(differences(:,3));
                if dvSortNew(1,1)<3
                    v1row = cat(1,v1row,dvSortOrig(1,1));
                end
                if i == 1 && size(v1row,1)>1       
                v1 = (rNDB(best,1:2) - rNDB(v1row(end,1)))/count;
                count = count+1;
                end
        end
    end
    
    clear v2row
    v2row = best;
    for i = 1:2
        count = 1;
        dv2 = rNDB(best,1:2);
        while (dv2(1,1) < rROIBounds(1,3) && dv2(1,1) > rROIBounds(1,1) && dv2(1,2) < rROIBounds(1,4) && dv2(1,2) > rROIBounds(1,2)) == 1

                clear differences
                dv2 = dv2 + v2*(-1^i);
                dv22(1:size(rNDB,1),1)=dv2(1,1);
                dv22(1:size(rNDB,1),2)=dv2(1,2);
                differences(:,1:2) = (rNDB(:,1:2) - dv22(:,1:2));
                differences(:,3) = sqrt(differences(:,1).^2+differences(:,2).^2);
                [dvSortNew,dvSortOrig] = sort(differences(:,3));
                if dvSortNew(1,1)<3
                    v2row = cat(1,v2row,dvSortOrig(1,1));
                end
                if i == 1 && size(v2row,1)>1       
                v2 = (rNDB(best,1:2) - rNDB(v2row(end,1)))/count;
                count = count+1;
                end

        end
    end
    v1row = unique(v1row);
    v2row = unique(v2row);
    
    [v1A,v1B] = fitLine3D(rNDB(v1row,1),rNDB(v1row,2),rNDB(v1row,3));
    [v2A,v2B] = fitLine3D(rNDB(v2row,1),rNDB(v2row,2),rNDB(v2row,3));
    
    for i=1:size(v1row,1)
        v1row(i,2) = norm(cross(v1B-v1A,rNDB(v1row(i,1),1:3)'-v1A))/norm(v1B-v1A);
    end
    
    for i=1:size(v2row,1)
        v2row(i,2) = norm(cross(v2B-v2A,rNDB(v2row(i,1),1:3)'-v2A))/norm(v2B-v2A);
    end
    v1mean = mean(v1row(:,2));
    v2mean = mean(v2row(:,2));
    %%
    if v1mean<v2mean
        rowV = (v1B-v1A)';
        better = 1;
    else
        rowV = (v2B-v2A)';
        better = 2;
    end
    
    d1 = (v1B-v1A)';
    d2 = (v2B-v2A)';
    %%
    v1B = v1B+(d1'*100);
    v1A = v1A-(d1'*100);
    v2B = v2B+(d2'*100);
    v2A = v2A-(d2'*100);
    
    viewDetections(rROI,raw,[0,0],image.TDots);
    hold on
    if better == 1
    plot3([v1A(1) v1B(1)]/raw.dataKey(9,1),[v1A(2) v1B(2)]/raw.dataKey(9,1),[v1A(3) v1B(3)],'g','LineWidth',20)
    plot3([v2A(1) v2B(1)]/raw.dataKey(9,1),[v2A(2) v2B(2)]/raw.dataKey(9,1),[v2A(3) v2B(3)],'r','LineWidth',20)
    else
    plot3([v1A(1) v1B(1)]/raw.dataKey(9,1),[v1A(2) v1B(2)]/raw.dataKey(9,1),[v1A(3) v1B(3)],'r','LineWidth',20)
    plot3([v2A(1) v2B(1)]/raw.dataKey(9,1),[v2A(2) v2B(2)]/raw.dataKey(9,1),[v2A(3) v2B(3)],'g','LineWidth',20)
    end
    hold off
    %% Display v1row and v2row
    
%     figure
%     hold on
%     scatter3(rNDB(:,1),rNDB(:,2),rNDB(:,3))
%     scatter3(rNDB(v1row(:,1),1),rNDB(v1row(:,1),2),rNDB(v1row(:,1),3))
%     scatter3(rNDB(v2row(:,1),1),rNDB(v2row(:,1),2),rNDB(v2row(:,1),3))
%     scatter3(0,0,0)
%     hold off
end