 close all; clear;
 set(0,'defaultfigurecolor',[1 1 1])
%Analyzing Trajectories from FIJI input or from Custom Code

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 1.) Inputs and Sorting the Raw Data From Excel File
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.1 Loading Data and Background Images for Overlays%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auto = questdlg('Attempt to Automate?',...
    'Automated inputs?','Yes','No','Yes');

if strcmp(auto,'Yes') == 1
    autoChk = 1;
else
    autoChk = 0;
end
disp('1.1 Loading Tracking Outputs')

[num,dataKey] = InputSelector(autoChk);

files = dir('*.tif'); %Check Directory for default filenames
filePath = strcat(cd,'\');
%Open Black Image
for k = 1:length(files)
    current=files(k).name;
    check(k)=strcmp(current(end-8:end),'black.tif');
end
loc=find(check);
if size(loc,1)==1
imageBlack = imread(files(loc(1)).name);
else
[nameBlackFile,filePath] = uigetfile('*.tif','Select a Black Image of the Correct Dimensions');
imageBlack = imread([filePath,nameBlackFile]);
end
 


if autoChk == 0
w = questdlg('Use a Transmitted image for overlays?',...
    'Transmitted Image (Optional)','Yes','No','Yes');
waitfor(w);
else
    w = 'Yes';
end
if strcmp(w,'Yes') == 1
    clear check
    for k = 1:length(files)
    current=files(k).name;
    if length(current)>=26 
    check(k)=strcmp(current(end-25:end),'Transmitted Cell Image.tif');
    end
    end
    loc=find(check);
    if size(loc,1)==1
    imageTrans = imread(files(loc(1)).name);
    else
    [nameTransFile,filePath] = uigetfile('*.tif','Select Transmitted Image for Overlay');
    imageTrans = imread([filePath,nameTransFile]);
    end
else
    imageTrans = imageBlack;
end

%Open Fluorescent Image
if autoChk == 0
w = questdlg('Use a Fluorescent image for overlays?',...
    'Fluorescent Image (Optional)','Yes','No','No');
waitfor(w);
else
    w = 'No';
end
if strcmp(w,'Yes') == 1
    clear check
    for k = 1:length(files)
    current=files(k).name;
    if length(current)>=26 
    check(k)=strcmp(current(end-26:end),'Fluorescent Cell Image.tif');
    end
    end
    loc=find(check);
    if size(loc,1)==1
    imageFluor = imread(files(loc(1)).name);
    else
    [nameFluorFile,filePath] = uigetfile('*.tif','Select Fluorescent Image for Overlay');
    imageFluor = imread([filePath,nameFluorFile]);
    end
else
    imageFluor = imageBlack;
end

%open Cell Area Binary Image
if autoChk == 0
w = questdlg('Input Binary Cell Outline?',...
    'Binary Outline','Yes','No','Yes');
else
    w = 'Yes';
end
if strcmp(w,'Yes') == 1
    
    clear check
    for k = 1:length(files)
    current=files(k).name;
    if length(current)>=15 
    check(k)=strcmp(current(end-14:end),'Binary Mask.tif');
    end
    end
    loc=find(check);
    if size(loc,1)==1
    imageArea= imread(files(loc(1)).name);
    imageArea = imresize(imageArea,(size(imageBlack,1)/size(imageArea,1)));
    else
    [nameAreaFile,filePath] = uigetfile('*.tif','Select a Thresholded Image of the Cell Area');
    imageArea = imread([filePath,nameAreaFile]);
    imageArea = imresize(imageArea,(size(imageBlack,1)/size(imageArea,1)));
    end
if size(imageArea,1) ~= size(imageBlack,1) || size(imageArea,2) ~= size(imageBlack,2)
    [nameAreaFile,filePath] = uigetfile('*.tif','Select a Thresholded Image of the Cell Area');
    imageArea = imread([filePath,nameAreaFile]);
end
else  
    imageArea = imageBlack==0;
end

%Open processed image stack of dots
    for k = 1:length(files)
    current=files(k).name;
    check(k)=strcmp(current(end-6:end),'roi.tif');
    end
    loc=find(check);
    if size(loc,1)==1
    roiStack = getImages(files(loc(1)).name);
    else
    roiStack = getImages();
    end

imageBorders = ones(size(imageArea,1),size(imageArea,2));
% bLimits = round(6/dataKey(9,1));
% imageBorders(end-bLimits:end,:) = 0;
% imageBorders(1:bLimits,:) = 0;
% imageBorders(:,end-bLimits:end) = 0;
% imageBorders(:,1:bLimits) = 0;


outputs = OutputSelector(autoChk);

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
clear book1
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
%8 = Global Tilt Corrected dX
%9 = Global Tilt Corrected dY
%10 = Global Tilt Corrected Magnitude
%11 = Mean local dx tilt
%12 = Mean local dy tilt
%13 = Local Tilt Corrected dx
%14 = Local Tilt Corrected dy
%15 = Local Tilt Corrected Magnitude
%16 = Index %15 is above cmCutoff (section 4.6)
%17 = Local Tilt Corrected dx above cutoff (using idx 16 as a mask)
%18 = Local Tilt Corrected dy above cutoff
%19 = Local Tilt Corrected Magnitude above cutoff

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
%11= Global Tilt Corrected Final dx
%12= Global Tilt Corrected Final dy
%13= Global Tilt Corrected Final Magnitude
%14= Predicted Top Surface
%15= Deviation from top Surface
%16= Alternative Predicted Top Surface
%17= Deviation from Alternative Top Surface
%18= Local Tilt Corrected Final dx
%19= Local Tilt Corrected Final dy
%20= Local Tilt Corrected Final Magnitude
%21= Book1 Index 16 is >0 for at least 3 consecutive frames

%HANDLING XYZ DATA
skipCount = 0;
missingNo = 0;
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
        book2((i-skipCount)-missingNo,3) = startFrame;
        book2((i-skipCount)-missingNo,4) = endFrame;
        % this fills in the upper portion of obj matrix with zeros if the first
        % frame is not 0
        
        book1(1,startFrame:endFrame,(i-skipCount)-missingNo) = tempObj(:,dataKey(1,1)).*dataKey(7,1);
        book1(2,startFrame:endFrame,(i-skipCount)-missingNo) = tempObj(:,dataKey(2,1)).*dataKey(7,1);
        book1(6,startFrame:endFrame,(i-skipCount)-missingNo) = tempObj(:,dataKey(5,1));
%      else
%          startFrame = 1;
%          endFrame = 1;
%          book2(i-skipCount,3) = startFrame;
%          book2(i-skipCount,4) = endFrame;
        end
    else
        missingNo = missingNo +1;
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
%% Find non-deformed pillars with no previous info
noCellTrajIni = 0;
SE = strel('disk',80);
imageAreaDil = imerode(imageArea,SE);
for i = 1:numTraj
    %if it is in black region
    if imageAreaDil(round(book2(i,2)),round(book2(i,1)))~=0
        if noCellTrajIni == 0
            noCellTrajIni = i;
        else
        noCellTrajIni = cat(1,noCellTrajIni,i);
        end
    end
end
imshow(imageAreaDil)

%% Tilt Correction
[noiseBook,noiseStats,sumIndFinal] = tiltCorrection(roiStack,imageTrans,book1,book2,noCellTrajIni);


%%
for i = 1:totalNumFrames
    book1(8,i,:) = ((book1(3,i,:)) - noiseBook(i,2)).*(book1(3,i,:)~=0);
    book1(9,i,:) = ((book1(4,i,:)) - noiseBook(i,3)).*(book1(4,i,:)~=0);
    book1(10,i,:) = ((book1(8,i,:).^2)+(book1(9,i,:).^2)).^.5;
    noiseBook(i,6) = mean(book1(8,i,sumIndFinal));
    noiseBook(i,7) = mean(book1(9,i,sumIndFinal));
    noiseBook(i,8) = mean(book1(10,i,sumIndFinal));
end

%%
bins = 0:.025:3;
errorHist = figure;
hold on
histogram(book1(10,:,sumIndFinal)*dataKey(9,1),bins,'normalization','probability')
histogram(book1(5,:,sumIndFinal)*dataKey(9,1),bins,'normalization','probability')
savefile = [filePath '\ErrorHistNoStress.tif'];
export_fig(errorHist,savefile);


errorHist2 = figure;
hold on
histogram((book1(10,:,:))*dataKey(9,1),bins,'normalization','probability')
histogram((book1(5,:,:))*dataKey(9,1),bins,'normalization','probability')
savefile = [filePath '\HistDisplacements.tif'];
export_fig(errorHist2,savefile);
%%
for i = 1:numTraj
    
    book2(i,5) = book1(3,book2(i,4),i);
    book2(i,6) = book1(4,book2(i,4),i);
    book2(i,7) = book1(1,book2(i,4),i);
    book2(i,8) = book1(2,book2(i,4),i);
    book2(i,9) = max(book1(5,:,i));%(book2(i,5).^2 + book2(i,6).^2).^0.5;
    book2(i,11) = book1(8,book2(i,4),i);
    book2(i,12) = book1(9,book2(i,4),i);
    book2(i,13) = (book2(i,11).^2 + book2(i,12).^2).^0.5;
end
%% Find all pillars in non-deformed regions
clear imageArea2
imageArea2 = zeros(size(imageArea,1)+100,size(imageArea,2)+100);
imageArea2(51:size(imageArea,1)+50,51:size(imageArea,2)+50) = imageArea;
for i = 1:numTraj
    if book2(i,13) > (1/dataKey(9,1))
        imageArea2((round(book2(i,2))):(round(book2(i,2))+100),(round(book2(i,1))):(round(book2(i,1))+100))=0;
    end
end
imageArea2 = imcrop(imageArea2,[51,51,size(imageArea,2)-1,size(imageArea,1)-1]);
noCellTraj = 0;
for i = 1:numTraj
    %if it is in black region
    if imageArea2(round(book2(i,2)),round(book2(i,1)))~=0
        if noCellTraj == 0
            noCellTraj = i;
        else
        noCellTraj = cat(1,noCellTraj,i);
        end
    end
end
%%
clear book4
%book4 contains the closest 500 pillars to index pillar(dim 1) that do not
%fall within the boundaries of the cell specified in var 'imageArea'

% find mean deviation of 500 closest in 'noCellTraj' variable
for i = 1:numTraj
    clear tempDistances
    if size(noCellTraj,1)>100
    tempDistances(1:size(noCellTraj,1),1) = ((book2(noCellTraj,1)-book2(i,1)).^2.+(book2(noCellTraj,2)-book2(i,2)).^2).^0.5;
    [tempDistances2,origOrder] = sortrows(tempDistances);
    book4(i,1:100) = noCellTraj(origOrder(1:100,1),1); %record the closest 500 pillars
    else
    tempDistances(1:size(noCellTraj,1),1) = ((book2(:,1)-book2(i,1)).^2.+(book2(:,2)-book2(i,2)).^2).^0.5;
    [~,origOrder] = sortrows(tempDistances);
    book4(i,1:size(noCellTraj,1)) = origOrder(1:size(noCellTraj,1),1); %record all pillars
    end
end
clear tempDistances
%%
clear book5
book5 = book1;
book5(book5(3:4,:,:) == 0) = NaN;

for j = 1:totalNumFrames
for i = 1:numTraj
    if book1(3,j,i) == 0 || book1(4,j,i) == 0
    book1(11,j,i) = 0;
    book1(12,j,i) = 0;
    book1(13,j,i) = 0;
    book1(14,j,i) = 0;
    book1(15,j,i) = 0;
    else
    book1(11,j,i) = nanmean(book5(3,j,book4(i,:)));
    book1(12,j,i) = nanmean(book5(4,j,book4(i,:)));
    book1(13,j,i) = book1(3,j,i)-book1(11,j,i);
    book1(14,j,i) = book1(4,j,i)-book1(12,j,i);
    book1(15,j,i) = ((book1(13,j,i).^2)+(book1(14,j,i).^2)).^.5;
    end
end
end

%%

for i = 1:numTraj

    book2(i,18) = book1(13,book2(i,4),i);
    book2(i,19) = book1(14,book2(i,4),i);
    book2(i,20) = max(book1(15,:,i));%(book2(i,18).^2 + book2(i,19).^2).^0.5;
  
end

% for i = 1:numTraj
%     if book1(13,book2(i,4)-1,i) == 0 || book1(14,book2(i,4)-1,i) == 0
%     book2(i,18) = max(book1(13,:,i));
%     book2(i,19) = max(book1(14,:,i));
%     book2(i,20) = (book2(i,18).^2 + book2(i,19).^2).^0.5;
%     else
%     book2(i,18) = book1(13,book2(i,4)-1,i);
%     book2(i,19) = book1(14,book2(i,4)-1,i);
%     book2(i,20) = (book2(i,18).^2 + book2(i,19).^2).^0.5;
%     end
% end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 3.) Identifying Neighborhood Trajectories for Each Trajectory
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% disp('3.1 Identifying Pillar Neighbor Groups')
% 
% clear tempDistances
% book3 = zeros(numTraj,50);
% cm3 = zeros(numTraj,50,totalNumFrames);
% cm4 = zeros(numTraj,totalNumFrames);
% interpXq = linspace(1,size(book1,2),size(book1,2)*3);
% maxDistance = (max(max(book1(5,:,:))));
% for i = 1:numTraj
%     tempDistances(1:numTraj,1) = ((book2(:,1)-book2(i,1)).^2.+(book2(:,2)-book2(i,2)).^2).^0.5;
%     tempDistances(1:numTraj,2) = linspace(1,numTraj,numTraj);
%     tempDistances = sortrows(tempDistances);
%     book3(i,1:50) = tempDistances(1:50,2); %record the closest 50 pillars
%     if book2(i,9) < (maxDistance) %filter by deformation limits here and in following if statement
%         cm4(i,1:totalNumFrames) = book1(6,1:totalNumFrames,i);
%     else
%         cm4(i,1:totalNumFrames) = NaN;
%     end
%     
%     for j = 1:50
%         if book2(book3(i,j),9) < (maxDistance) % *following if statement*
%             cm3(i,j,1:totalNumFrames) = book1(6,1:totalNumFrames,book3(i,j)); %record intensity value if pillar has low deformation
%         else
%             cm3(i,j,1:totalNumFrames) = NaN;
%         end
%     end
% end
% 
% cm4(cm4==0) = NaN;
% totalAverage = zeros(1,totalNumFrames);
% for i = 1:totalNumFrames
%     totalAverage(1,i) = mean(cm4(:,i),'omitnan');
% end
% totalAverageInterp(1,:) = interp1(linspace(1,size(book1,2),size(book1,2)),conv(totalAverage(1,:),gausswin(6),'same'),interpXq,'spline');
% 

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
disp('4.2 Creating Color Map for Vector Plot Image Overlays')
[cm1,cm2,cmD,cmDS,colorMap,colorScheme] = createColorMap(book1,book2,outputs);




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.3 Using Intensity Values to Extract Z-Information%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: Calculates average distance between PSFs in a single pillar.
%Outputs stored as variable "yDataDiffAverage"

%*****Not really very useful as of 2017*****

if ismember(10,outputs) == 1
    disp('4.3 Fitting Marker Plane')
    planeFit(book1,book2,totalNumFrames,filePath)
end




%%
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

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4.6 Data Filters and Cutoffs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.6.1 Persistence - The goal is to remove shear deformations that "occur"
% in less than 3 consecutive frames. Whether a shear deformation "occurs"
% will depend on whether the displacement is greater than some threshold.
disp('4.6 Thresholding Data Using a Displacement Magnitude-Based Cutoff')
cmCutoff = 3; %The first colormap group to be used in filtered data sets
book1(16,:,:) = book1(15,:,:)>2;
cmConsec = ones(3,1);
for i = 1:size(book1,3)
    if sum(book1(16,:,i),2)<3
        book2(i,21) = 0;
    else
        book2(i,21) = 1;
    end
    
    if book2(i,21) == 1
        temp = book1(16,:,i);
        for j = 2:size(temp,2)-1
            temp2(j,1) = temp(1,j)+temp(1,j-1)+temp(1,j+1);
        end
        if max(temp2)<3
            book2(i,21) = 0;
        end
    end
end
%%
filterSet = find(book2(:,21)>0);
for i = 1:size(filterSet,1);
    book1(17:19,:,filterSet(i,1)) = book1(13:15,:,filterSet(i,1));
end
if size(filterSet,1)<1
    book1(17:19,:,:) = 0;
end

filterMask = book1(19,:,:)>2;
book1(17,:,:) = book1(17,:,:).*filterMask;
book1(18,:,:) = book1(18,:,:).*filterMask;
book1(19,:,:) = book1(19,:,:).*filterMask;
book2(:,22) = book2(:,21).*book2(:,20);
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%5.)IMAGE OUTPUTS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
disp('5.0 Creating Centroid and Vector Plot Overlay Image Outputs')
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
            quiver(book1(1,f,cm1(cm1(:,i,f)>0,i,f)),book1(2,f,cm1(cm1(:,i,f)>0,i,f)),book1(13,f,cm1(cm1(:,i,f)>0,i,f)),book1(14,f,cm1(cm1(:,i,f)>0,i,f)),...
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
            quiver(book2(cm2(cm2(:,i)>0,i),1),book2(cm2(cm2(:,i)>0,i),2),book2(cm2(cm2(:,i)>0,i),18),book2(cm2(cm2(:,i)>0,i),19),0,'color',[colorMap(i,1:3)]);
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
        
        quiver(book2(cm2(cm2(:,i)>0,i),1),book2(cm2(cm2(:,i)>0,i),2),book2(cm2(cm2(:,i)>0,i),18),book2(cm2(cm2(:,i)>0,i),19),0,'color',[colorMap(i,1:3)]);%(i^2)/((cmD/1.5)^2) also (i^1.3)/((cmD/2)^1.3) ...(i^2)/((cmD/1.5)^2)
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
        quiver(book2(cm2(cm2(:,i)>0,i),1),book2(cm2(cm2(:,i)>0,i),2),book2(cm2(cm2(:,i)>0,i),18),book2(cm2(cm2(:,i)>0,i),19),0,'color',[colorMap(i,1:3)]);%(i^1.3)/((cmD/2)^1.3)
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

%% ***WIP*** Create FE Mesh

if ismember(11,outputs) == 1
    disp('6.0 Creating FE Mesh')
   meshbook(); 
end
%% Create Heat Map of Shear Deformations
if ismember(15,outputs) == 1
[imageHeatXY,vqXY,vqXYtotal,imageHeatNaN,imageHeatXYColor]=customHeatMap(book1,book2,imageBlack,dataKey,outputs,filePath);
end
SE = strel('disk',30);
vqXYBinary = imdilate(vqXY ==0,SE);
BinaryFile = [filePath,'HeatMaps\Single\','Binary_Shear.tif'];
        imwrite(vqXYBinary,BinaryFile);
%%
imageHeatXY2 = imresize(vqXY,[size(imageArea,1) size(imageArea,2)]);
imageHeatXY2(isnan(imageHeatXY2)) = 0;
%%
clear  imageHeatXYTotalScale imageHeatXYtotal2
for i = 1:size(vqXYtotal,3)
imageHeatXYtotal2(:,:,i) = imresize(vqXYtotal(:,:,i),[size(imageArea,1) size(imageArea,2)]);
end
imageHeatXYtotal2(isnan(imageHeatXYtotal2)) = 0;

%%
[fitTiltPlaneMicrons,fitTiltPlanePixels] = tiltPlane(noiseBook,dataKey,imageArea);

%%
clear topSurface fitSurface
[topSurface ,fitSurface] = findSurface(book1,book2,cm2,cmCutoff,imageArea,imageBorders,dataKey);
save('fitSurface.mat','fitSurface')
%%
figure
%plot the filter pillars at the top frame that they reach
scatter3(0,0,0)
hold on
scatter3(topSurface(:,3),topSurface(:,2),topSurface(:,4));
%plot(fitSurface{1})
plot(fitSurface{2})
 xlim([0 size(roiStack,1)])
 ylim([0 size(roiStack,2)])
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
book2(i,16) = feval(fitSurface{2},book1(2,book2(i,4),i),book1(1,book2(i,4),i));
book2(i,15) = 0;
book2(i,17) = 0;
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
book2(cellSurface(i,1),17) = book2(cellSurface(i,1),4)-book2(cellSurface(i,1),16);
end


interpSurface{1} = fit([(book2(:,2) + book2(:,12)),(book2(:,1) + book2(:,11))],(book2(:,14)+book2(:,15)),'lowess','Span',.01);

%% Calculate Z-Displacement at Surface (Quick Method 2)
for i = 1:numTraj
book2(i,14) = feval(fitSurface{1},book1(2,book2(i,4),i),book1(1,book2(i,4),i));
book2(i,16) = feval(fitSurface{2},book1(2,book2(i,4),i),book1(1,book2(i,4),i));
book2(i,15) = book2(i,4)-book2(i,14);
book2(i,17) = book2(i,4)-book2(i,16);
end


interpSurface{1} = fit([(book2(:,2) + book2(:,12)),(book2(:,1) + book2(:,11))],(book2(:,14)+book2(:,15)),'lowess','Span',.005);

%%
% close all
figure
quiver3(book2(:,2),book2(:,1),book2(:,14),book2(:,12),book2(:,11),book2(:,15))
hold on
plot3(0,0,0)
plot(interpSurface{1})
xlim([0 size(roiStack,1)])
ylim([0 size(roiStack,2)])
zlim([0 size(roiStack,3)]) 
hold off

%%
if ismember(16,outputs) == 1
    %%
[imageHeatZ,vqZ] = customHeatMapZ(book2,imageBlack,dataKey,outputs,filePath);
end
%% Isolating Cell-Body Normal Forces
%Filter Z-Deformation to accept only normal deformation within the cell's boundary
imageHeatZScale = size(imageArea,1)/size(vqZ,1);
imageHeatZ2 = imresize(vqZ,[size(imageArea,1) size(imageArea,2)]);
imageHeatZ3 = double(imageHeatZ2) .* double(imageArea==0);
imageHeatZ3(isnan(imageHeatZ3)) = 0;


%% Relate Cell Spread Area to Displacements
clear imageAreaProps imageArea3
book1(find(isnan(book1))) = 0;
imageArea3(:,:) = logical(imageArea==0);
imageAreaProps = regionprops(imageArea3,'centroid','Area','Eccentricity','Perimeter','MajorAxisLength','MinorAxisLength');
propsArea = sum(sum(imageArea==0));
propsSumDisp = sum(sum(book1(19,:,:)));
propsMeanDisp = mean(mean(book1(19,:,:)));
if propsArea > 0
propsCircularity = ((sum(cat(1,imageAreaProps.Perimeter)))^2 )/(4*(pi*(sum(cat(1,imageAreaProps.Area)))));
ratioArea = propsSumDisp/propsArea;
ratioPerimeter = propsSumDisp/cat(1,imageAreaProps.Perimeter);
else
    clear imageAreaProps
    imageAreaProps = struct('centroid',0,'Area',0,'Eccentricity',0,'Perimeter',0,'MajorAxisLength',0,'MinorAxisLength',0);
%     imageAreaProps.Perimeter =0;
%     imageAreaProps.Area = 0;
%     imageAreaProps(1,1).MajorAxisLength = 0;
%     imageAreaProps(1,1).MinorAxisLength = 0;
%     imageAreaProps(1,1).Eccentricity = 0;
    propsCircularity = 0;
    ratioArea = 0;
    ratioPerimeter = 0;
end


areaTxt = fopen('Area-Disp Relationship.txt','wt');
fprintf(areaTxt,strcat('Code Date:07/13/2017','\n'));

fprintf(areaTxt,strcat( num2str(sum(cat(1,imageAreaProps.Area))) , '\n'));
fprintf(areaTxt,strcat( num2str(propsSumDisp) , '\n'));
fprintf(areaTxt,strcat( num2str(propsMeanDisp) ,'\n'));
fprintf(areaTxt,strcat( num2str(sum(cat(1,imageAreaProps.Perimeter))) ,'\n'));
fprintf(areaTxt,strcat( num2str(propsCircularity) ,'\n'));
fprintf(areaTxt,strcat( num2str(imageAreaProps(1,1).Eccentricity) ,'\n'));
fprintf(areaTxt,strcat(num2str(imageAreaProps(1,1).MajorAxisLength) , '\n'));
fprintf(areaTxt,strcat(num2str(imageAreaProps(1,1).MinorAxisLength) , '\n'));
fprintf(areaTxt,strcat( num2str(ratioArea) ,'\n'));
fprintf(areaTxt,strcat(num2str(ratioPerimeter) ,'\n'));
fprintf(areaTxt,strcat(num2str(sum(sum(imageHeatXY2))) , '\n'));
fprintf(areaTxt,strcat(num2str(sum(sum(sum(imageHeatXYtotal2)))) , '\n'));
fprintf(areaTxt,strcat(num2str(sum(sum(imageHeatZ3))) , '\n'));
fprintf(areaTxt,strcat(num2str(max(max(imageHeatXY2))) , '\n'));
fprintf(areaTxt,strcat('\n'));
fprintf(areaTxt,strcat('Above is Without Text for Copy Paste','\n'));
fprintf(areaTxt,strcat('\n'));
fprintf(areaTxt,strcat( num2str(sum(cat(1,imageAreaProps.Area))) ,',    Cell Area in Sq. Pixels: ', '\n'));
fprintf(areaTxt,strcat( num2str(propsSumDisp) ,',   Sum of Pillar Displacements in Linear Pixels: ', '\n'));
fprintf(areaTxt,strcat( num2str(propsMeanDisp) ,',  Average Displacement in Linear Pixels: ', '\n'));
fprintf(areaTxt,strcat( num2str(sum(cat(1,imageAreaProps.Perimeter))) ,',   Perimeter: ', '\n'));
fprintf(areaTxt,strcat( num2str(propsCircularity) ,',   Circularity: ', '\n'));
fprintf(areaTxt,strcat( num2str(imageAreaProps(1,1).Eccentricity) ,',   Eccentricity', '\n'));
fprintf(areaTxt,strcat(num2str(imageAreaProps(1,1).MajorAxisLength) ,',   Major Axis of Ellipse', '\n'));
fprintf(areaTxt,strcat(num2str(imageAreaProps(1,1).MinorAxisLength) ,',   Minor Axis of Ellipse', '\n'));
fprintf(areaTxt,strcat( num2str(ratioArea) ,',  Area Ratio: ', '\n'));
fprintf(areaTxt,strcat(num2str(ratioPerimeter) ,',  Perimeter Ratio: ', '\n'));
fprintf(areaTxt,strcat(num2str(sum(sum(imageHeatXY2))) ,',  Sum of Interpolated XY Displacements in Linear Microns (Surface): ', '\n'));
fprintf(areaTxt,strcat(num2str(sum(sum(sum(imageHeatXYtotal2)))) ,',    Sum of Interpolated XY Displacements in Linear Microns (Stack): ', '\n'));
fprintf(areaTxt,strcat(num2str(sum(sum(imageHeatZ3))) ,',   Sum of Interpolated Normal Displacements in Linear Microns (Cell Boundary): ', '\n'));
fprintf(areaTxt,strcat(num2str(max(max(imageHeatXY2))) ,',   Max of Interpolated XY Displacements in Linear Microns (Surface): ', '\n'));

fclose(areaTxt);
%% Create folder for profile views and save relevant data
filePath = cd;
folderName = 'Profile Data';
mkdir(filePath,folderName)
save('Profile Data\vqXY.mat','vqXY')
save('Profile Data\HeatMapXY.mat','imageHeatXYColor')

 %% Calculate Z-Displacement at Surface (Quick Method 2 in microns)

% %interpSurface{1} = fit([(book2(:,2) + book2(:,12))*dataKey(9,1),(book2(:,1) + book2(:,11))*dataKey(9,1)],(book2(:,14)+book2(:,15))*0.4,'lowess','Span',.01);
% interpSurface{1} = fit([(book2(:,2) + book2(:,12))*dataKey(9,1),(book2(:,1) + book2(:,11))*dataKey(9,1)],(book2(:,16)+book2(:,17))*0.4,'lowess','Span',.01);
% 
% close all
% figure
% quiver3(book2(:,2)*dataKey(9,1),book2(:,1)*dataKey(9,1),book2(:,14)*0.4,book2(:,12)*dataKey(9,1),book2(:,11)*dataKey(9,1),book2(:,15)*0.4,'r')
% hold on
% plot3(0,0,0)
% plot(interpSurface{1})
% xlim([0 size(roiStack,1)]*dataKey(9,1))
% ylim([0 size(roiStack,2)]*dataKey(9,1))
% zlim([0 size(roiStack,3)]*0.4)
% hold off
%%
% 
% figure
% hold on
% quiver3(book2(:,2),book2(:,1),book2(:,14),book2(:,12),book2(:,11),book2(:,15))
% plot(interpSurface{2})
% plot3(0,0,0)

% %% Calculate Area and Displacement
% for j = 1:5
% 
% [nameAreaFile,filePath] = uigetfile('*.tif','Select a Thresholded Image of the Cell Area');
% imageArea2 = imread([filePath,nameAreaFile]);
% Area2 = sum(sum(imageArea2==0));
% 
% [nameCOIFile,filePath] = uigetfile('*.tif','Select a Thresholded Image of the Cell Forces');
% imageForceCOI = imread([filePath,nameCOIFile]);
% 
% 
% forces = 0;
% for i = 1:numTraj
%     %if it is in black region
%     if imageForceCOI(round(book2(i,2)),round(book2(i,1)))==0
%         if forces == 0
%             forces = i;
%         else
%         forces = cat(1,forces,i);
%         end
%     end
% end
% 
% for i = 1:totalNumFrames
%     forcesZ(i,1) = sum(book1(10,i,forces));
% end
% 
% sumDisp = sum(forcesZ);
% averageDisp = sumDisp/Area2;
% 
% cd(filePath)
% 
% areaTxt = fopen('Area-Disp Relationship.txt','wt');
% fprintf(areaTxt,strcat('Cell Area in Sq. Pixels: ', num2str(Area2) , '\n'));
% fprintf(areaTxt,strcat('Sum of Displacements in Linear Pixels: ', num2str(sumDisp) , '\n'));
% fprintf(areaTxt,strcat('Ratio: ', num2str(averageDisp) , '\n'));
% 
% for i = 1:totalNumFrames
%     fprintf(areaTxt,strcat('Sum of Displacements in Linear Pixels for Frame ',num2str(i),': ', num2str(forcesZ(i,1)) , '\n'));
% end
% 
% fclose(areaTxt);
% 
% AreaDispCOIs(1,j) = Area2;
% AreaDispCOIs(2,j) = sumDisp;
% AreaDispCOIs(3,j) = averageDisp;
% AreaDispCOIs(4:(3+totalNumFrames),j) = forcesZ;
% 
% end
%
% txtFileName = 'AreaDispPerCell.txt';
%     if exist(txtFileName,'file')
%         delete(txtFileName)
%     end
%     dlmwrite(txtFileName,AreaDispCOIs);
% 

%%
% close all
% i = 1000;
% plot(1:1:totalNumFrames,book1(15,:,i))
% hold on
% plot(1:1:totalNumFrames,book1(11,:,i))
% plot(1:1:totalNumFrames,book1(12,:,i))
% figure
% hold on
% plot(1:1:totalNumFrames,book1(13,:,i)) 
% plot(1:1:totalNumFrames,book1(14,:,i))
% figure
% plot(1:1:totalNumFrames,book1(3,:,i))
% hold on
% plot(1:1:totalNumFrames,book1(4,:,i))
disp('Trajectories Program has Completed Successfully')