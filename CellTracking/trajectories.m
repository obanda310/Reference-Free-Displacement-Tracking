close all; clear;
%Analyzing Trajectories from FIJI input or from Custom Code

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 1.) Inputs and Sorting the Raw Data From Excel File
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.1 Loading Data and Background Images for Overlays%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1.1 Loading Tracking Outputs')

[num,dataKey] = InputSelector();

[nameBlackFile,filePath] = uigetfile('*.tif','Select a Black Image of the Correct Dimensions');
imageBlack = imread([filePath,nameBlackFile]);

w = questdlg('Use a Fluorescent image for overlays?',...
    'Fluorescent Image (Optional)','Yes','No','No');
waitfor(w);
if strcmp(w,'Yes') == 1
    [nameFluorFile,filePath] = uigetfile('*.tif','Select Fluorescent Image for Overlay');
    imageFluor = imread([filePath,nameFluorFile]);
else
    imageFluor = imageBlack;
end

w = questdlg('Use a Transmitted image for overlays?',...
    'Fluorescent Image (Optional)','Yes','No','No');
waitfor(w);
if strcmp(w,'Yes') == 1
    [nameTransFile,filePath] = uigetfile('*.tif','Select Transmitted Image for Overlay');
    imageTrans = imread([filePath,nameTransFile]);
else
    imageTrans = imageBlack;
end

w = questdlg('Input Binary Cell Outline?',...
    'Binary Outline','Yes','No','No');
if strcmp(w,'Yes') == 1
    [nameAreaFile,filePath] = uigetfile('*.tif','Select a Thresholded Image of the Cell Area');
    imageArea = imread([filePath,nameAreaFile]);
    imageArea = imresize(imageArea,(size(imageBlack,1)/size(imageArea,1)));
else
    imageArea = imageBlack==0;
end

roiStack = getImages();

imageBorders = ones(size(imageArea,1),size(imageArea,2));
% bLimits = round(6/dataKey(9,1));
% imageBorders(end-bLimits:end,:) = 0;
% imageBorders(1:bLimits,:) = 0;
% imageBorders(:,end-bLimits:end) = 0;
% imageBorders(:,1:bLimits) = 0;


outputs = OutputSelector();

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.2 Misc Input Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1.2 Loading Variables')
numIndices = 7; %for Section 1 and 2 (number of elements in book1)
numTraj = max(num(:,dataKey(4,1))); %Number of Trajectories
totalNumFrames = max(num(:,dataKey(3,1)))+dataKey(8,1); %Maximum number of frames observable for any one object

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 2.) Building book1 and book2 from the raw data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

disp('2.1 Creating book1 and book2')

% book1 = zeros(numIndices,totalNumFrames,numTraj);
% book2 = zeros(numTraj,10);

%Notes: The majority of indices are added in later sections, but listed
%here for reference.

%Index List for Book1 - Frame Dependent Values
%1 = Raw X Centroid Location of Traj in Specified Frame
%2 = Raw Y Centroid Location of Traj in Specified Frame
%3 = Raw dX Centroid Location of Traj in Specified Frame from First Frame
%4 = Raw dY Centroid Location of Traj in Specified Frame from First Frame
%5 = Magnitude of displacement
%6 = Intensity in current frame (column 14 from TrackMate output)
%7 = Gauss filtered intensity values (created/used in section 4)
%8 = Tilt Corrected dX
%9 = Tilt Corrected dY
%10 = Tilt Corrected Magnitude

%Index List for Book2 - Frame Independent Values
%1 = Raw X Centroid Location of Traj in First Frame
%2 = Raw Y Centroid Location of Traj in First Frame
%3 = First Frame that a Traj Appears
%4 = Last Frame that a Traj Appears
%5 = Value of book1 index 3 (see above) in Last Frame that Traj Appears
%6 = Value of book1 index 4 (see above) in Last Frame that Traj Appears
%7 = Value of book1 index 1 (see above) in Last Frame that Traj Appears
%8 = Value of book1 index 2 (see above) in Last Frame that Traj Appears
%9 = Maximum Magnitude of displacement
%10= Numeric ID of Pillar
%11= Tilt Corrected Final dx
%12= Tilt Corrected Final dy
%13= Tilt Corrected Final Magnitude
%14= Predicted Top Surface
%15= Deviation from top Surface

%HANDLING XYZ DATA
skipCount = 0;
for i = 1:numTraj
    % Here we build a book of pages (3D array) with the data for a single
    % object/trajectory per page. Most of the code is to ensure that each
    % page is the same size matrix as the next.
    
    % stores data to 'tempObj' pertaining to all frames of one
    % object/trajectory
    tempObj = num(num(:,dataKey(4,1))==i,:);
    % number of frames that the current object appears in.
    numFrames = size(tempObj,1);
    if numFrames > 0
        % the first frame that an object appears in (tracking software starts
        % at 0
        startFrame = min(tempObj(:,dataKey(3,1)))+dataKey(8,1);
        % the last frame that an object appears in
        endFrame = (startFrame+numFrames-1);
        
        if imageBorders(round(tempObj(1,3)*dataKey(7,1)),round(tempObj(1,2)*dataKey(7,1)))==0
            skip = 1;
            skipCount = skipCount+1;
        else
            skip = 0;
        end
        
        if skip == 0
        book2(i-skipCount,3) = startFrame;
        book2(i-skipCount,4) = endFrame;
        % this fills in the upper portion of obj matrix with zeros if the first
        % frame is not 0
        
        book1(1,startFrame:endFrame,i-skipCount) = tempObj(:,dataKey(1,1)).*dataKey(7,1);
        book1(2,startFrame:endFrame,i-skipCount) = tempObj(:,dataKey(2,1)).*dataKey(7,1);
        book1(6,startFrame:endFrame,i-skipCount) = tempObj(:,dataKey(5,1));
%     else
%         startFrame = 1;
%         endFrame = 1;
%         book2(i-skipCount,3) = startFrame;
%         book2(i-skipCount,4) = endFrame;
        end
    end
    if i == 1 || i == 5000 || i == 10000 || i == 15000 || i==20000 || i==25000
        disp(['Progress: ' num2str(i) ' of ' num2str(numTraj)])
    end
end
numTraj = size(book1,3);

disp('2.2 Storing XY Displacements')
for i = 1:numTraj
    book2(i,1) = book1(1,book2(i,3),i); %initial x
    book2(i,2) = book1(2,book2(i,3),i); %initial y
    book1(3,book2(i,3):book2(i,4),i) = book1(1,book2(i,3):book2(i,4),i) - book2(i,1); %frame specific dx
    book1(4,book2(i,3):book2(i,4),i) = book1(2,book2(i,3):book2(i,4),i) - book2(i,2); %frame specific dy
    book1(5,:,i) = (book1(3,:,i).^2 + book1(4,:,i).^2).^0.5; %frame specific total displacement
    book2(i,10) = i; %pillar ID
end

% Tilt Correction
[noiseBook,sumIndFinal] = tiltCorrection(roiStack,imageTrans,book1,book2);

for i = 1:totalNumFrames
    book1(8,i,:) = ((book1(3,i,:)) - noiseBook(i,2)).*(book1(3,i,:)~=0);
    book1(9,i,:) = ((book1(4,i,:)) - noiseBook(i,3)).*(book1(4,i,:)~=0);
    book1(10,i,:) = ((book1(8,i,:).^2)+(book1(9,i,:).^2)).^.5;
    noiseBook(i,6) = mean(book1(8,i,sumIndFinal));
    noiseBook(i,7) = mean(book1(9,i,sumIndFinal));
    noiseBook(i,8) = mean(book1(10,i,sumIndFinal));
end

errorHist = figure;
hold on
histogram(book1(10,:,sumIndFinal),'normalization','probability')
histogram(book1(5,:,sumIndFinal),'normalization','probability')
savefile = [filePath '\ErrorHistNoStress.tif'];
export_fig(errorHist,savefile);


errorHist2 = figure;
hold on
histogram(book1(10,:,:),'normalization','probability')
histogram(book1(5,:,:),'normalization','probability')
savefile = [filePath '\HistDisplacements.tif'];
export_fig(errorHist2,savefile);

for i = 1:numTraj
    
    book2(i,5) = book1(3,book2(i,4),i);
    book2(i,6) = book1(4,book2(i,4),i);
    book2(i,7) = book1(1,book2(i,4),i);
    book2(i,8) = book1(2,book2(i,4),i);
    book2(i,9) = (book2(i,5).^2 + book2(i,6).^2).^0.5;
    book2(i,11) = book1(8,book2(i,4),i);
    book2(i,12) = book1(9,book2(i,4),i);
    book2(i,13) = (book2(i,5).^2 + book2(i,6).^2).^0.5;
end


%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 3.) Identifying Neighborhood Trajectories for Each Trajectory
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

disp('3.1 Identifying Pillar Neighbor Groups')

clear tempDistances
book3 = zeros(numTraj,50);
cm3 = zeros(numTraj,50,totalNumFrames);
cm4 = zeros(numTraj,totalNumFrames);
interpXq = linspace(1,size(book1,2),size(book1,2)*3);
maxDistance = (max(max(book1(5,:,:))));
for i = 1:numTraj
    tempDistances(1:numTraj,1) = ((book2(:,1)-book2(i,1)).^2.+(book2(:,2)-book2(i,2)).^2).^0.5;
    tempDistances(1:numTraj,2) = linspace(1,numTraj,numTraj);
    tempDistances = sortrows(tempDistances);
    book3(i,1:50) = tempDistances(1:50,2); %record the closest 50 pillars
    if book2(i,9) < (maxDistance) %filter by deformation limits here and in following if statement
        cm4(i,1:totalNumFrames) = book1(6,1:totalNumFrames,i);
    else
        cm4(i,1:totalNumFrames) = NaN;
    end
    
    for j = 1:50
        if book2(book3(i,j),9) < (maxDistance) % *following if statement*
            cm3(i,j,1:totalNumFrames) = book1(6,1:totalNumFrames,book3(i,j)); %record intensity value if pillar has low deformation
        else
            cm3(i,j,1:totalNumFrames) = NaN;
        end
    end
end

cm4(cm4==0) = NaN;
totalAverage = zeros(1,totalNumFrames);
for i = 1:totalNumFrames
    totalAverage(1,i) = mean(cm4(:,i),'omitnan');
end
totalAverageInterp(1,:) = interp1(linspace(1,size(book1,2),size(book1,2)),conv(totalAverage(1,:),gausswin(6),'same'),interpXq,'spline');


%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 4.)Identifying Deviations From Zero-State in Later Frames
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.1 Storing the Maximum Displacement Values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: Some trajectories do not persist all the way to the top frame, and
%so they would not be visible in the outputs in section 5 unless their
%maximum value is used for those overlays.
disp('4.1 Storing Maximum Displacement Values')
%moved to previous section



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.2 Creating a Color Map for Quiver Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cm1,cm2,cmD,cmDS,colorMap,colorScheme] = createColorMap(book1,book2,outputs);




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.3 Using Intensity Values to Extract Z-Information%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: Calculates average distance between PSFs in a single pillar.
%Outputs stored as variable "yDataDiffAverage"

%*****Not really very useful as of 2017*****

if ismember(10,outputs) == 1
    planeFit(book1,book2,totalNumFrames,filePath)
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.4 Using Intensity Values to Extract Z-Information%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For use later
zScale = .4; %microns
highestObjectInZ = (size(book1,2)-1)* zScale;


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
    pillarPlot(book1,book2,book3,cm3,roiStack,totalNumFrames,interpBook,degreeInterp);
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
    
    finalLoc(:,1:2,:) = finalLoc(:,1:2,:)*dataKey(9,1);
    
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
            if finalLocSorted(i,1,j)<finalLocSorted(rowStart,1,j)+10*dataKey(9,1)
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



%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%5.)IMAGE OUTPUTS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%5.2 Drawing Zero-State Displacement Fields on Black Background%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if ismember(2,outputs) == 1
    folderName = strcat('Black Image_',colorScheme,' Quiver Overlays');
    mkdir(filePath,folderName)
    for f = 1:totalNumFrames        % number of z-slices
        blackOverlay = figure('Position',[0 0 1000 1000]);
        imshow(imageBlack,[])
        hold on
        
        for i = 1:cmD
            quiver(book1(1,f,cm1(cm1(:,i,f)>0,i,f)),book1(2,f,cm1(cm1(:,i,f)>0,i,f)),book1(3,f,cm1(cm1(:,i,f)>0,i,f)),book1(4,f,cm1(cm1(:,i,f)>0,i,f)),...
                0,'color',[colorMap(i,1:3)]);
            hold on
        end
        hold off
        savefile = [filePath '\' folderName '\Black Background Overlay' colorScheme ' ' num2str(f) '.tif'];
        if ismember(6,outputs) == 1
            export_fig(blackOverlay,savefile,'-native');
        else
            export_fig(blackOverlay,savefile);
        end
        close
    end
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%5.3 Plotting Centroids on a Black Background%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: This section plots centroids (as in 5.3) on a black background for
%processing into heat maps in imageJ to view 3D information in a 2D image.

if ismember(4,outputs) == 1
    folderName = strcat( 'Centroid_',colorScheme,' Black Image Overlays');
    mkdir(filePath, folderName)
    for f = 1:totalNumFrames
        centroidsOnly = figure('Position',[0 0 1000 1000]);
        imshow(imageBlack,[])
        hold on
        for i = 1:cmD
            
            xTemp = squeeze(book1(1,f,cm1(cm1(:,i,f)>0,i,f)));
            
            yTemp = squeeze(book1(2,f,cm1(cm1(:,i,f)>0,i,f)));
            
            plot(xTemp,yTemp,'.','MarkerSize',30,'Color',[colorMap(i,1:3)]);
            
            hold on
        end
        hold on
        hold off
        savefile = [filePath '\' folderName '\Centroids on Frame ' colorScheme ' ' num2str(f) '.tif'];
        if ismember(6,outputs) == 1
            export_fig(centroidsOnly,savefile,'-native');
        else
            export_fig(centroidsOnly,savefile);
        end
        close
    end
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%5.4 DEBUGGING:Plotting Corrected X,Y Coordinate Fields%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: This sections is really helpful for debugging. It shows x0,y0 positions
% as well as current X,Y coordinates for each trajectory on each frame. It
% also plots the final displacement quiver field.

if ismember(3,outputs) == 1
    debugImage = figure('Position',[0 0 1000 1000]);
    imshow(imageTrans,[])
    hold on
    xTemp = squeeze(book1(1,1,:));
    yTemp = squeeze(book1(2,1,:));
    plot(xTemp,yTemp,'r.','MarkerSize',10);
    hold on
    if ismember(7,outputs) == 0
        for f = 1:totalNumFrames
            xTemp = squeeze(book1(1,f,:));
            yTemp = squeeze(book1(2,f,:));
            plot(xTemp,yTemp,'b.','MarkerSize',3);
            hold on
        end
    end
    hold on
    if ismember(8,outputs) == 0
        for i = 1:cmD
            quiver(book2(cm2(cm2(:,i)>0,i),1),book2(cm2(cm2(:,i)>0,i),2),book2(cm2(cm2(:,i)>0,i),11),book2(cm2(cm2(:,i)>0,i),12),0,'color',[colorMap(i,1:3)]);
            hold on
        end
        
    end
    hold off
    savefile = [filePath '\Debug Image ' colorScheme '.tif'];
    if ismember(6,outputs) == 1
        export_fig(debugImage,savefile,'-native');
    else
        export_fig(debugImage,savefile);
    end
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%5.5 Plotting Transmitted Quiver Overlays v1%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: 5.5 Saves files, 5.6 does not.

if ismember(5,outputs) == 1
    trajOverlay = figure('Position',[0 0 1000 1000]);
    imshow(imageFluor,[])
    hold on
    for i = 1:cmD
        
        quiver(book2(cm2(cm2(:,i)>0,i),1),book2(cm2(cm2(:,i)>0,i),2),book2(cm2(cm2(:,i)>0,i),11),book2(cm2(cm2(:,i)>0,i),12),0,'color',[colorMap(i,1:3)]);%(i^2)/((cmD/1.5)^2) also (i^1.3)/((cmD/2)^1.3) ...(i^2)/((cmD/1.5)^2)
        hold on
    end
    %quiver(book1(10,1,:),book1(11,1,:),book1(12,totalNumFrames,:),book1(13,totalNumFrames,:),0,'g');
    hold off
    savefile = [filePath '\Fluorescent Overlay ' colorScheme '.tif'];
    if ismember(6,outputs) == 1
        export_fig(trajOverlay,savefile,'-native');
    else
        export_fig(trajOverlay,savefile);
    end
end

%%
if ismember(5,outputs) == 1
    transmittedOverlay = figure('Position',[0 0 1000 1000]);
    imshow(imageTrans,[])
    hold on
    for i = 1:cmD
        quiver(book2(cm2(cm2(:,i)>0,i),1),book2(cm2(cm2(:,i)>0,i),2),book2(cm2(cm2(:,i)>0,i),11),book2(cm2(cm2(:,i)>0,i),12),0,'color',[colorMap(i,1:3)]);%(i^1.3)/((cmD/2)^1.3)
        hold on
    end
    %quiver(book1(10,1,:),book1(11,1,:),book1(12,totalNumFrames,:),book1(13,totalNumFrames,:),0,'g');
    hold off
    savefile = [filePath '\Transmitted Overlay ' colorScheme '.tif'];
    if ismember(6,outputs) == 1
        export_fig(transmittedOverlay,savefile,'-native');
    else
        export_fig(transmittedOverlay,savefile);
    end
end

%%
if ismember(11,outputs) == 1
    
    % Create Mesh from 2D shear data
    
    clear MeshBook MeshList
    clear book4
    
    book4 = book1;
    book4(~book4) = NaN;
    
    MeshBook = book4(1:2,1:25,1:121);
    MeshBook(2,:,:) = MeshBook(2,:,:)*-1;
    node = 0;
    
    for i = 1:size(MeshBook,2)
        for j = 1:size(MeshBook,3)
            node = node + 1;
            if isnan(MeshBook(1,i,j)) == 1
                MeshBook(1,i,j) = MeshBook(1,i-1,j);
            end
            if isnan(MeshBook(2,i,j)) == 1
                MeshBook(2,i,j) = MeshBook(2,i-1,j);
            end
            MeshBook(3,i,j) = i;
            MeshBook(4,i,j) = node;
            
            
            MeshList(2,node) = MeshBook(1,i,j);
            MeshList(3,node) = MeshBook(2,i,j);
            MeshList(4,node) = MeshBook(3,i,j);
            MeshList(1,node) = MeshBook(4,i,j);
            
            MeshList2(2,node) = MeshBook(1,1,j);
            MeshList2(3,node) = MeshBook(2,1,j);
            MeshList2(4,node) = MeshBook(3,i,j);
            MeshList2(1,node) = MeshBook(4,i,j);
        end
    end
    
    %Building Elements
    %Start with the first frame
    clear EleList EleList2
    Element = 0;
    for i = 1:size(MeshBook,3)
        clear current diff
        current(1:2,:) = MeshBook(1:2,1,:);
        diff(1,:) = current(1,:) - MeshBook(1,1,i);
        diff(2,:) = current(2,:) - MeshBook(2,1,i);
        diff(3,:) = sqrt(diff(1,:).^2.+diff(2,:).^2);
        diff(diff<-5) = NaN;
        diff(diff>22) = NaN;
        for j = 1:size(diff,2)
            if isnan(diff(1,j))==1 || isnan(diff(2,j))== 1 || isnan(diff(3,j)) ==1 || diff(3,j) == 0
                diff(1:3,j) = 0;
            end
        end
        
        clear current2
        current2 = find(diff(3,:));
        
        if size(current2,2) == 3
            current3 = diff(1:3,current2);
            mesh3 = current2(find(current3(3,:)==max(current3(3,:))));
            mesh4 = current2(find(current3(2,:)==min(current3(2,:))));
            mesh2 = current2(find(current3(1,:)==min(current3(1,:))));
            for k = 1:(size(MeshBook,2)-1)
                Element = Element + 1;
                
                EleList(i,1,k) = MeshBook(4,k,i); %5
                EleList(i,2,k) = MeshBook(4,k,mesh2); %6
                EleList(i,3,k) = MeshBook(4,k,mesh3); %7
                EleList(i,4,k) = MeshBook(4,k,mesh4); %8
                EleList(i,5,k) = MeshBook(4,k+1,i); %1
                EleList(i,6,k) = MeshBook(4,k+1,mesh2); %2
                EleList(i,7,k) = MeshBook(4,k+1,mesh3); %3
                EleList(i,8,k) = MeshBook(4,k+1,mesh4); %4
                
                EleList(i,9,k) = Element;
                EleList2(2:9,Element) = EleList(i,1:8,k);
                EleList2(1,Element) = Element;
            end
        end
        
        
    end
    
    %Write the deformed mesh
    
    meshTxt = fopen('mesh.txt','wt');
    nodesFormat = '<node id=" %d "> %f, %f, %f </node>\n';
    fprintf(meshTxt,'<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_spec version="2.5">\n<Geometry>\n<Nodes name="Part1">\n');
    fprintf(meshTxt,nodesFormat,MeshList(1:4,:));
    fprintf(meshTxt,'</Nodes>\n<Elements type="hex8" mat="1" name="part1">\n');
    elementsFormat = '<elem id=" %d "> %d, %d, %d, %d, %d, %d, %d, %d </elem>\n';
    fprintf(meshTxt,elementsFormat,EleList2(1:9,:));
    fprintf(meshTxt,'</Elements>\n</Geometry></febio_spec>');
    fclose(meshTxt);
    
    %Write the undeformed Mesh
    
    meshTxt = fopen('meshUD.txt','wt');
    nodesFormat = '<node id=" %d "> %f, %f, %f </node>\n';
    fprintf(meshTxt,'<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_spec version="2.5">\n<Geometry>\n<Nodes name="Part1">\n');
    fprintf(meshTxt,nodesFormat,MeshList2(1:4,:));
    fprintf(meshTxt,'</Nodes>\n<Elements type="hex8" mat="1" name="part1">\n');
    elementsFormat = '<elem id=" %d "> %d, %d, %d, %d, %d, %d, %d, %d </elem>\n';
    fprintf(meshTxt,elementsFormat,EleList2(1:9,:));
    fprintf(meshTxt,'</Elements>\n</Geometry></febio_spec>');
    fclose(meshTxt);
    
    %Display the nodes from the deformed mesh
    
    figure
    hold on
    for i = 1:size(MeshBook,2)
        scatter3(MeshBook(1,i,:),MeshBook(2,i,:),i*ones(1,size(MeshBook,3)),'.','g');
        scatter3(MeshBook(1,1,:),MeshBook(2,1,:),i*ones(1,size(MeshBook,3)),'.','b');
    end
    hold off
    
end


%%
customHeatMap(book1,book2,imageBlack,dataKey,outputs,filePath)

%%
cmArea = sum(sum(imageArea==0));
cmCutoff = ceil(2*(mean(mean(book1(10,:,sumIndFinal))))/cmDS);
for i = 1:totalNumFrames
    for j = cmCutoff:cmD
        temp(i,j) = sum(book1(10,i,cm1(cm1(:,j,i)>0,j,i)));
    end
end
sumDisp = sum(sum(temp));
averageDisp = sumDisp/cmArea;

areaTxt = fopen('Area-Disp Relationship.txt','wt');
fprintf(areaTxt,strcat('Cell Area in Sq. Pixels: ', num2str(cmArea) , '\n'));
fprintf(areaTxt,strcat('Sum of Displacements in Linear Pixels: ', num2str(sumDisp) , '\n'));
fprintf(areaTxt,strcat('Ratio: ', num2str(averageDisp) , '\n'));
fclose(areaTxt);


%%
[fitTiltPlaneMicrons,fitTiltPlanePixels] = tiltPlane(noiseBook,dataKey,imageArea);

%%
clear topSurface fitSurface
[topSurface ,fitSurface] = findSurface(book1,book2,cm2,cmCutoff,imageArea,imageBorders);

%%
figure
%plot the filter pillars at the top frame that they reach
scatter3(0,0,0)
hold on
scatter3(topSurface(:,3),topSurface(:,2),topSurface(:,4));
plot(fitSurface{1})
%plot(fitTiltPlanePixels)
hold off
%%
figure
imshow(imageArea)
hold on
scatter(topSurface(:,2),topSurface(:,3),5,'r')
hold off

%% Calculate Z-Displacement at Surface (Quick Method)
clear cellSurface
for i = 1:numTraj
book2(i,14) = feval(fitSurface{1},book1(2,book2(i,4),i),book1(1,book2(i,4),i));
book2(i,15) = 0;
end

cellSurface = zeros(1,1);
for i = 1:numTraj
    %if it is under the cell
    if imageArea(round(book2(i,2)),round(book2(i,1)))==0 && imageArea(round(book2(i,8)),round(book2(i,7)))==0       
        cellSurface = cat(1,cellSurface,i);
    end
end
%shift cells up 1 to get rid of initial zero
cellSurface(1,:) = [];

for i = 1:size(cellSurface,1)
book2(cellSurface(i,1),15) = book2(cellSurface(i,1),4)-book2(cellSurface(i,1),14);
end


%interpSurface{1} = fit([(book2(:,2) + book2(:,12)),(book2(:,1) + book2(:,11))],(book2(:,14)+book2(:,15)),'lowess','Span',.01);

%% Calculate Z-Displacement at Surface (Quick Method 2)
for i = 1:numTraj
book2(i,14) = feval(fitSurface{1},book1(2,book2(i,4),i),book1(1,book2(i,4),i));
book2(i,15) = book2(i,4)-book2(i,14);
end


%interpSurface{1} = fit([(book2(:,2) + book2(:,12)),(book2(:,1) + book2(:,11))],(book2(:,14)+book2(:,15)),'lowess','Span',.005);

%%
% close all
% figure
% quiver3(book2(:,2),book2(:,1),book2(:,14),book2(:,12),book2(:,11),book2(:,15))
% hold on
% plot3(0,0,0)
% plot(interpSurface{1})
% xlim([0 size(roiStack,1)])
% ylim([0 size(roiStack,2)])
% zlim([0 size(roiStack,3)])
% hold off

%%

%% Calculate Z-Displacement at Surface (Quick Method 2 in microns)

interpSurface{1} = fit([(book2(:,2) + book2(:,12))*dataKey(9,1),(book2(:,1) + book2(:,11))*dataKey(9,1)],(book2(:,14)+book2(:,15))*0.4,'lowess','Span',.005);

close all
figure
quiver3(book2(:,2)*dataKey(9,1),book2(:,1)*dataKey(9,1),book2(:,14)*0.4,book2(:,12)*dataKey(9,1),book2(:,11)*dataKey(9,1),book2(:,15)*0.4,'r')
hold on
plot3(0,0,0)
plot(interpSurface{1})
xlim([0 size(roiStack,1)]*dataKey(9,1))
ylim([0 size(roiStack,2)]*dataKey(9,1))
zlim([0 size(roiStack,3)]*0.4)
hold off
%%
% 
% figure
% hold on
% quiver3(book2(:,2),book2(:,1),book2(:,14),book2(:,12),book2(:,11),book2(:,15))
% plot(interpSurface{2})
% plot3(0,0,0)

%%
disp('Trajectories Program has Completed Successfully')