%% Old code. Z positioning is now found using 3D centroid detections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.4 Using Intensity Values to Extract Z-Information%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For use later
zScale = .4; %microns
highestObjectInZ = (size(book1,2)-1)* zScale;

if ismember(17,outputs) == 1
disp('4.4 Extracting Z-Information Based on through-Z Intensity Profiles')
%Interpolation to increase resolution of intensity data points
clear interpBook interpX
degreeInterp = 12;
interpBook = zeros(size(book1,3),size(book1,2)*degreeInterp);
interpX = linspace(1,size(book1,2),size(book1,2));
interpXq = linspace(1,size(book1,2),size(book1,2)*degreeInterp);
for i = 1:numTraj
    book1(7,:,i) = conv(book1(6,:,i),gausswin(6),'same');
    interpBook(i,:) = interp1(interpX,book1(7,:,i),interpXq,'spline');
end

%Find the z-dimension peaks based on intensity
clear zPeaks
zPeaks = zeros(numTraj,10);
for i = 1:numTraj
    clear pks loc truePeaks
    [pks,loc] = findpeaks(interpBook(i,:));
    truePeaks = find(pks>max(pks)*.25);
    zPeaks(i,1:size(loc(truePeaks),2)) = loc(truePeaks);
    
    %     if zPeaks(i,1) <7 % remove false peaks created at bottom limit of image stack
    %         zPeaks(i,1:8) = zPeaks(i,2:9);
    %     end
    
    zPeaks(i,10) = nnz(zPeaks(i,1:9));
end

%add a psuedo first peak for data alignment
% for i = 1:size(zPeaks,1)
%     if  zPeaks(i,1)>degreeInterp*5
%         zPeaks(i,2:9) = zPeaks(i,1:8);
%         zPeaks(i,1) = 1;
%     end
% end


%"Statistically" find problem pillars with incorrect or missing peaks
zPeakAvgNo = floor(mean(zPeaks(:,10)));
clear zPeaksSpacing zSpacingIssues
zPeaksSpacing(1:size(zPeaks,1),1) = 30; %Space filler larger than issue threshold below
zPeaksSpacing(:,2:zPeakAvgNo+1) = zPeaks(:,3:zPeakAvgNo+2)-zPeaks(:,2:zPeakAvgNo+1);
for i = 2:size(zPeaksSpacing,2)
    for j = 1:size(zPeaks,1)
        if abs(zPeaksSpacing(j,i))<degreeInterp*7
            zPeaks(j,i) = round(mean(zPeaks(j,i:i+1)));
            zPeaks(j,(i+1):8)=zPeaks(j,(i+2):9);
        end
    end
end

%Correct shifts in data due to uneven recognition of objects
zPeaksAverage = mean(zPeaks(:,2:zPeakAvgNo),2);
% for i = 1:size(zPeaks,1)
%     if

%     end
% end


zPeaksSpacing = zPeaks(:,2:zPeakAvgNo+2)-zPeaks(:,1:zPeakAvgNo+1);
zPeaksSpacing(zPeaksSpacing<0) = NaN;
zPeaksSpacing(zPeaksSpacing==0) = NaN;
zPeaksAvgSpacing = mean(mean(zPeaksSpacing,'omitnan'),'omitnan');
zPeaksStDSpacing = std(mean(zPeaksSpacing,'omitnan'),'omitnan');
%find pillars with spacing outside of mean+/- 10*STD


%Create a list of substitute peaks for pillars with problems
clear zSubPeaks zAvgClose10 zAvgClose10Interp
zAvgClose50 = zeros(size(cm3,1),size(cm3,3));
zAvgClose50Interp = zeros(size(cm3,1),size(cm3,3)*degreeInterp);
zAvgClose50(:,:) = mean(cm3(:,1:50,:),2,'omitnan');
zSubPeaks = zeros(numTraj,10);
for i = 1:numTraj
    zAvgClose50Interp(i,:) = interp1(linspace(1,size(book1,2),size(book1,2)),conv(zAvgClose50(i,:),gausswin(6),'same'),interpXq,'spline');
    clear pks loc truePeaks
    [pks,loc] = findpeaks(zAvgClose50Interp(i,:));
    truePeaks = find(pks>max(pks)/4);
    
    zSubPeaks(i,1:size(loc(truePeaks),2)) = loc(truePeaks);
    
    %     if zPeaks(i,1) <7 % remove false peaks created at bottom limit of image stack
    %         zPeaks(i,1:8) = zPeaks(i,2:9);
    %     end
end

[zPeaksIssues,~] = find(zPeaks(:,1:zPeakAvgNo)==0);
%Replace Problem Pillars 'zPeaksIssues' with an average pillar
for i = 1:size(zPeaksIssues,1)
    zPeaks(zPeaksIssues(i,1),1:zPeakAvgNo) = zSubPeaks(zPeaksIssues(i,1),1:zPeakAvgNo);
end


zPeaksDifferences = zPeaks(2:end,1:zPeakAvgNo)-zPeaks(1:end-1,1:zPeakAvgNo);
[a,~] = find(abs(zPeaksDifferences)>((2/zScale)*degreeInterp));
zPeaksIssues2 = find(histc(a,1:numTraj)>1);

%Replace Problem Pillars 'zPeaksIssues' with an average pillar
for i = 1:size(zPeaksIssues2,1)
    zPeaks(zPeaksIssues2(i,1),1:zPeakAvgNo) = zSubPeaks(zPeaksIssues2(i,1),1:zPeakAvgNo);
end


%Convert zPeaks to a frame value based on previous interpolation
zPeaks(:,1:9) = zPeaks(:,1:9)./degreeInterp;

%Create a layer based on the maximum
for i = 1:numTraj
    zPeaks(i,zPeakAvgNo+1) = book2(i,4);
end
zPeakAvgNo = zPeakAvgNo +1;



%%
if ismember(9,outputs) == 1
    pillarPlot(book1,book2,book3,cm3,image.ROIstack,totalNumFrames,interpBook,degreeInterp);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.5 Estimating XY Coordinates of Interpolated Z Positions%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if zPeakAvgNo > 2
    clear zPeaksFrames zPeaksTopWeights zPeaksBottomWeights finalLoc
    
    zPeaksFrames = floor(zPeaks);
    zPeaksTopWeights = zPeaks-zPeaksFrames;
    zPeaksBottomWeights = 1-zPeaksTopWeights;
    
    finalLoc = zeros(numTraj,4,zPeakAvgNo);
    
    for i = 1:numTraj
        for j = 1:zPeakAvgNo-1
            %Structure is (Top Frame Weight Fraction * Top Frame) + (Bottom Frame Weight
            %Fraction * Bottom Frame) All Divided by the number of parts (the degree of
            %interpolation stored in 'degreeInterp'). The whole thing is just a
            %weighted average.
            
            %first check if top frame and bottom frame are non-zero
            %for xPos
            if book1(1,zPeaksFrames(i,j)+1,i)>0 && zPeaksFrames(i,j)>0
                finalLoc(i,1,j) =  (((zPeaksTopWeights(i,j)*degreeInterp)*book1(1,zPeaksFrames(i,j)+1,i))+((zPeaksBottomWeights(i,j)*degreeInterp)*book1(1,zPeaksFrames(i,j),i)))/degreeInterp; %x position
            else %if both are not nonzero use last available position
                finalLoc(i,1,j) = book1(1,book2(i,4),i);
            end
            %for yPos
            if book1(2,zPeaksFrames(i,j)+1,i)>0 && zPeaksFrames(i,j)>0
                finalLoc(i,2,j) =  (((zPeaksTopWeights(i,j)*degreeInterp)*book1(2,zPeaksFrames(i,j)+1,i))+((zPeaksBottomWeights(i,j)*degreeInterp)*book1(2,zPeaksFrames(i,j),i)))/degreeInterp;%y position
            else
                finalLoc(i,2,j) = book1(2,book2(i,4),i);
            end
            
            finalLoc(i,3,j) =  zPeaks(i,j)*zScale;%z position
            finalLoc(i,4,j) = i;
        end
        
    end
    
    for i = 1:numTraj
        finalLoc(i,1,zPeakAvgNo) = book1(1,book2(i,4),i);
        finalLoc(i,2,zPeakAvgNo) = book1(2,book2(i,4),i);
        finalLoc(i,3,zPeakAvgNo) = zPeaks(i,zPeakAvgNo)*zScale;
        finalLoc(i,4,zPeakAvgNo) = i;
    end
    
    finalLoc(:,1:2,:) = finalLoc(:,1:2,:)*raw.dataKey(9,1);
    
    %Create ZX profile for Ryan
    clear finalLocSorted values order endFitData
    for i = 1:zPeakAvgNo
        [values, order] = sort(finalLoc(:,1,i));
        finalLocSorted(:,:,i) = finalLoc(order,:,i);
    end
    
    clear Xs Zs
    
    for j = 1:size(finalLocSorted,3)
        count1 = 1;
        count2 = 0;
        rowStart = 1;
        for i = 1:size(finalLocSorted,1)
            if finalLocSorted(i,1,j)<finalLocSorted(rowStart,1,j)+10*raw.dataKey(9,1)
                count2 = count2+1;
                Xs(count1,count2,j) = finalLocSorted(i,1,j);
                Zs(count1,count2,j) = finalLocSorted(i,3,j);
            else
                count1 = count1+1;
                count2 = 1;
                rowStart = i;
                Xs(count1,count2,j) = finalLocSorted(i,1,j);
                Zs(count1,count2,j) = finalLocSorted(i,3,j);
            end
        end
    end
    
    Xs(Xs==0) = NaN;
    Zs(Zs==0) = NaN;
    
    Xmeans = mean(Xs,2,'omitnan');
    Zmeans = mean(Zs,2,'omitnan');
    
    fitFraction =.1;
    leftFit = round(size(Xmeans,1)*fitFraction);
    rightFit = size(Xmeans,1) - leftFit;
    
    
    endsFitData = zeros(leftFit*2,2,zPeakAvgNo);
    for i = 1:zPeakAvgNo
        endsFitData(1:leftFit,1,i) = Xmeans(1:leftFit,1,i);
        endsFitData(1:leftFit,2,i) = Zmeans(1:leftFit,1,i);
        endsFitData(leftFit+1:(leftFit*2),1,i) = Xmeans(rightFit+1:end,1,i);
        endsFitData(leftFit+1:(leftFit*2),2,i) = Zmeans(rightFit+1:end,1,i);
    end
    endsFitDataTrunc = endsFitData(1:end-5,:,:);
    
    
    for i = 1:zPeakAvgNo
        fitObj{i} = polyfit(endsFitDataTrunc(:,1,i),endsFitDataTrunc(:,2,i),1);
    end
    
    %y = polyval(p,x) basic structure
    clear finalZPos
    for i = 1:zPeakAvgNo
        finalZPos(:,1,i) = (((max(endsFitDataTrunc(:,2,zPeakAvgNo))-polyval(fitObj{zPeakAvgNo},Xmeans(:,1,zPeakAvgNo))) + Zmeans(:,1,i)))-highestObjectInZ;
    end
    
    
    avgXZName = 'avgXZ.txt';
    dlmwrite(avgXZName,cat(2,Xmeans,Zmeans));
    
    
    %%
    if ismember(14,outputs) == 1
        %%
        %Show an XYZ representation of the dot positions
        figure
        hold on
        for i = 1:size(finalLoc,3)-1
            fitobject4{i} = fit([finalLoc(:,2,i),finalLoc(:,1,i)],finalLoc(:,3,i),'lowess','Span',0.05);
            scatter3(finalLoc(:,2,i),finalLoc(:,1,i),finalLoc(:,3,i),10);
            plot(fitobject4{i})
        end
        hold off
        %%
        %Show an XZ representation of the dot positions
        figure
        for i = 1:size(finalLoc,3)
            
            scatter(finalLoc(:,1,i),finalLoc(:,3,i));
            hold on
        end
        hold off
        
        %Show a flattened XZ representation of the average dot positions
        zReference = max(endsFitDataTrunc(:,2,zPeakAvgNo));
        figure
        for i = 1:size(Xmeans,3)
            plot(Xmeans(:,1,i),finalZPos(:,1,i));
            hold on
            plot(Xmeans(:,1,i),(polyval(fitObj{i},Xmeans(:,1,i)))+(zReference - polyval(fitObj{zPeakAvgNo},Xmeans(:,1,zPeakAvgNo)))-highestObjectInZ)
        end
        
        hold off
    end
    
end

end