%current scale 3.0769  (x4=12.3076)
close all; clear;
%Analyzing Trajectories from FIJI input

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 1.) Inputs and Sorting the Raw Data From Excel File
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.1 Loading Data and Background Images for Overlays%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[name,path] = uigetfile('*.xlsx','Select .xlsx File From Particle Tracker Output');
file = [path,name];
[num,txt,raw] = xlsread(file);

inputVar = InputSelector();
if strcmp(inputVar,'Mosaic') == 1
    xCol = 4;
    yCol = 5;
    fCol = 3;
    tCol = 2;
    totalCol = 12;
    pixelScale = 1;
elseif strcmp(inputVar,'TrackMate') == 1
    xCol = 6;
    yCol = 7;
    fCol = 10; %Frame
    tCol = 4; %Pillar/Trajectory
    intCol = 14; %Intensity information
    totalCol = 22; %Number of columns in spreadsheet
    startVar = 1;
    prompt = 'How many pixels per micron? Enter a decimal and press enter: ';
    pixelScale = input(prompt);
elseif strcmp(inputVar,'Custom Code') == 1
    xCol = 2;
    yCol = 3;
    fCol = 4;
    tCol = 5;
    intCol = 6;
    totalCol = 6;
    prompt = 'What was the scale factor print out of tracking.m? Check the command window. Enter a decimal and press enter: ';
    pixelScale = input(prompt);
    startVar = 0;
end
    
[trajFile,trajPath] = uigetfile('*.tif','Select Fluorescent Image for Overlay');
c = imread([trajPath,trajFile]);

[cellFile,cellPath] = uigetfile('*.tif','Select Transmitted Image for Overlay');
d = imread([cellPath,cellFile]);

[blackFile,blackPath] = uigetfile('*.tif','Select a Black Image of the Correct Dimensions');
e = imread([blackPath,blackFile]);

outputs = OutputSelector();

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.2 Input Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
attemptRowAdjust = 0; %Try to adjust for wobble from patterning process? 1 for yes
[resY,resX] = size(c); %for Section 5
scaleOutputVectors = 1; %for Section 5.1
trackErrorThreshold = .2; %for Section 4
numIndices = 30; %for Section 1 and 2 (number of elements in book1)
numTraj = max(num(:,tCol)); %Number of Trajectories
totalNumFrames = max(num(:,fCol))+startVar; %Maximum number of frames observable for any one object
global book1 
book1 = zeros(numIndices,totalNumFrames,numTraj); %Creates book1

%%

for i = 1:numTraj
    % Here we build a book of pages (3D array) with the data for a single
    % object/trajectory per page. Most of the code is to ensure that each
    % page is the same size matrix as the next.
    
    % stores data to 'tempObj' pertaining to all frames of one
    % object/trajectory
    tempObj = num(num(:,tCol)==i,:);
    % number of frames that the current object appears in.
    numFrames = size(tempObj,1);
    if numFrames > 0
    % the first frame that an object appears in (tracking software starts 
    % at 0
    startFrame = min(tempObj(:,fCol))+startVar;
    % the last frame that an object appears in
    endFrame = (startFrame+numFrames-1);
    book1(11,:,i) = startFrame;
    book1(12,:,i) = endFrame;
    % this fills in the upper portion of obj matrix with zeros if the first
    % frame is not 0

    book1(1,startFrame:endFrame,i) = tempObj(:,xCol).*pixelScale;
    book1(2,startFrame:endFrame,i) = tempObj(:,yCol).*pixelScale;
    book1(30,startFrame:endFrame,i) = tempObj(:,intCol);   
    else
        startFrame = 1;
        endFrame = 1;
        book1(11,:,i) = startFrame;
        book1(12,:,i) = endFrame;
    end
    if i == 1 || i == 5000 || i == 10000 || i == 15000 || i==20000 || i==25000
        progress = 'Section 1' 
    end
end

 %Here we use the temporary data stored in 'obj' to build 2 matrices to
 %describe the x and y positions in frame 1 for each object, as well as the
 %displacement from those matrices occurring in successive frames.



%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 2.) Building book1 from the raw data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%Notes: The majority of indices are added in later sections, but listed 
%here for reference.

%Index List for Book1
%1 = Raw X Centroid Location of Traj in Specified Frame
%2 = Raw Y Centroid Location of Traj in Specified Frame
%3 = Raw X Centroid Location of Traj in First Frame
%4 = Raw Y Centroid Location of Traj in First Frame
%5 = Raw dX Centroid Location of Traj in Specified Frame from First Frame
%6 = Raw dY Centroid Location of Traj in Specified Frame from First Frame
%7 = Rounded dX Centroid Location of Traj in Specified Frame from First Frame
%8 = Rounded dY Centroid Location of Traj in Specified Frame from First Frame
%9 = Row that current Traj is a member of
%10 = Nearest Right neighbor Traj to current Traj
%11 = First Frame that a Traj Appears
%12 = Last Frame that a Traj Appears
%13 = Mode of Rounded X Displacements of all Traj in Current Traj's Row
%14 = Mode of Rounded X Displacements of all Traj in Current Traj's Row
%15 = Mode Corrected X_0 Coordinate of Traj
%16 = Mode Corrected Y_0 Coordinate of Traj
%17 = Mode Corrected X Displacement of Traj in Current Frame from First 
%18 = Mode Corrected Y Displacement of Traj in Current Frame from First
%19 = Magnitude of displacement
%20 = Nearest Left neighbor Traj to current Traj
%21 = Mode Corrected X Coordinate of Traj in Current Frame
%22 = Mode Corrected Y Coordinate of Traj in Current Frame
%23 = Value of book1 index 17 (see above) in Last Frame that Traj Appears
%24 = Value of book1 index 18 (see above) in Last Frame that Traj Appears
%25 = Value of book1 index 21 (see above) in Last Frame that Traj Appears
%26 = Value of book1 index 22 (see above) in Last Frame that Traj Appears
%27 = Full page of x_0
%28 = Full page of y_0
%29 = Maximum Magnitude of displacement
%30 = Intensity in current frame (column 14 from TrackMate output)
for i = 1:numTraj

        book1(3,book1(11,1,i):book1(12,1,i),i) = book1(1,book1(11,1,i),i);
        book1(4,book1(11,1,i):book1(12,1,i),i) = book1(2,book1(11,1,i),i);
        book1(5,:,i) = book1(1,:,i) - book1(3,:,i);
        book1(6,:,i) = book1(2,:,i) - book1(4,:,i);
        book1(7,:,i) = round(book1(5,:,i),1);
        book1(8,:,i) = round(book1(6,:,i),1);
        book1(19,:,i) = (book1(5,:,i).^2 + book1(6,:,i).^2).^0.5;
        book1(27,:,i) = book1(3,book1(11,1,i),i);
        book1(28,:,i) = book1(4,book1(11,1,i),i);

        
end
%%


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 3.) Assigning Objects to Rows by Based on Neighbor Particles
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%CURRENTLY DISABLED

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%3.1 Identifying Nearest Left and Right Neighbors%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if attemptRowAdjust == 1
objAll = permute(squeeze([book1(3,1,:),resY-book1(4,1,:)]),[2 1]);
rows = zeros(2,1);
row = 1;
col = 1;
for i = 1:numTraj
    if book1(11,1,i) ~= 1
        i=i+1;
    else    
        thisObj = squeeze([book1(3,1,i),resY-book1(4,1,i)]);
        objCount = (1:numTraj)';
        unweightedDistances = bsxfun(@minus,objAll,thisObj);
        xDistances = unweightedDistances(:,1);
        yDistances = unweightedDistances(:,2);
        distances = sqrt((xDistances.^2+yDistances.^2));
        objData = [objCount,distances];
        sortDist = sortrows(objData,2);
        neighbors(1,1:4) = sortDist(2:5,1);
        nDistance = sortDist(2:5,2);
        for j = 1:4
            n(j,1:2) = objAll(neighbors(1,j),:);
            n(j,3) = neighbors(j);
        end
        a = bsxfun(@minus,n(:,1:2),thisObj);
        b = mean([min(abs(a)),max(abs(a))]);
        goodNeighbors = n(abs(a(:,2)) < b & a(:,1) > 0,:);
        if size(goodNeighbors,1) > 1
            sortN = sortrows(goodNeighbors,1);
            rightNeighbor = sortN(1,:);
        else
            rightNeighbor = goodNeighbors;
        end
        goodNeighborsForLeft = n(abs(a(:,2)) < b & a(:,1) < 0,:);
        if size(goodNeighborsForLeft,1) > 1
            sortLeftN = sortrows(goodNeighborsForLeft,1);
            leftNeighbor = sortLeftN(1,:);
        else
            leftNeighbor = goodNeighborsForLeft;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%3.2 Row Assignment Method 1: Using Distance%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %A threshold y distance is chosen. Any trajectories outside of the
        %threshold distance in y are considered part of a different row than
        %the current row.
        yThresholdFactor = 0.5;
        if i>1 && abs(yDistances(i-1,1)) < abs(yThresholdFactor*nDistance(1,1))
            col = col + 1;
            % Add one to the width of the 'rows' array to have a slot 
            % which will necessarily be empty, stored as 'emptyColumn'
            emptyColumn = size(rows,2)+1;
            % If there is an empty value to the left of 'emptyColumn' 
            % in 'trow' then,
                  while rows(book1(9,1,i-1),emptyColumn-1) == 0
                            emptyColumn = emptyColumn-1;
                  end
            % place the trajectory number into a new empty column
            rows(book1(9,1,i-1),emptyColumn) = i;
            book1(9,:,i) = row;
            % if a right neighbor doesn't exist for the current trajectory (i.e. 
            % if it is on the edge of the frame)
        elseif i>1 && abs(yDistances(i-1,1)) > abs(0.5*nDistance(1,1))
            row = row + 1;
            col = 1;
            rows(row,col) = i;
            book1(9,:,i) = row;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3.3 Row Assignment Method 2: Using Neighbor Trajectories%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        elseif isempty(rightNeighbor) == 1 && isempty(leftNeighbor) == 0   
            % adjust the current column to the next slot(to the right)
            col = col + 1;
            % place the current trajectory in that slot
            rows(row,col) = i;
            book1(9,:,i) = row;
            % then start a new row for the next trajectory in the list
            % (because we have reached an edge)
            row = row + 1;
            % start at the beginning of the new row for the next trajectory
            col = 1;
            book1(10,:,i)= 0;
            % if a left neighbor doesn't exist for the current trajectory 
            % (i.e. the first member of a row)
         elseif isempty(leftNeighbor)
            col = 1;    % adjust the current column to the first column
            rows(row,col) = i;  % place the current trajectory in that slot
            book1(9,:,i) = row;
            % start at the beginning of the new row for the next trajectory
            col = 1;
            book1(20,:,i)= 0;
         else
            % if we are NOT at the beginning or end of the row, we assume 
            % first that we MAY NOT be in the correct row, we check the 
            % current trajectory's right/left neighbor's row to see if they
            % exist 
            % is the right/left neighbor to current object already assigned a 
            % row? store as 'rowCheck'
            rowCheckRight = ismember(rightNeighbor(3),rows);
            rowCheckLeft = ismember(leftNeighbor(3),rows);
            if rowCheckRight == 1	% if right exists then
                % what row does the right neighbor belong to? Store as trow
                [trow,tcol] = find(rows == rightNeighbor(3));
                % Add one to the width of the 'rows' array to have a slot 
                % which will necessarily be empty, stored as 'emptyColumn'
                emptyColumn = size(rows,2)+1;
                % If there is an empty value to the left of 'emptyColumn' 
                % in 'trow' then,
                    while rows(trow,emptyColumn-1) == 0
                        emptyColumn = emptyColumn-1;
                    end
                % place the trajectory number into a new empty column
                rows(trow,emptyColumn) = i;
                book1(10,:,i)= rightNeighbor(3);
                book1(9,:,i) = trow;
                col = col + 1;
            elseif rowCheckLeft == 1	% if left exists then
                % what row does the left neighbor belong to? Store as trow
                [trow,tcol] = find(rows == leftNeighbor(3));
                % Add one to the width of the 'rows' array to have a slot 
                % which will necessarily be empty, stored as 'emptyColumn'
                emptyColumn = size(rows,2)+1;
                % If there is an empty value to the left of 'emptyColumn' 
                % in 'trow' then,
                    while rows(trow,emptyColumn-1) == 0
                        emptyColumn = emptyColumn-1;
                    end             
                % place the trajectory number into a new empty column
                rows(trow,emptyColumn) = i;
                book1(9,:,i) = trow;
                col = col + 1;             
            else
                % if the right neighbor is not on the list of rows yet, 
                % then move one column right and...
                col = col + 1;
                % ... add the current trajectory to the empty slot
                rows(row,col) = i;
                % add rightNeighbor to index book1 (index#10)

                book1(9,:,i) = row;
            end
            book1(10,:,i)= rightNeighbor(3);
            book1(20,:,i)= leftNeighbor(3);
        end
        if isempty(rightNeighbor)==0
        book1(10,:,i)= rightNeighbor(3);
        end
        if isempty(leftNeighbor)==0
        book1(20,:,i)= leftNeighbor(3);
        end
    end
end
end                                             
%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 4.)Identifying Deviations From Zero-State in Later Frames
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%CURRENTLY DISABLED

% Now that we have a list of rows of trajectories we need to find the
% error displacement for members of a particular row caused by 
% irregularities in the patterning process. 
%
% We are measuring from the zero frame to the current frame. 
%
% We can do this by rounding the x and y displacement of each value to the 
% nearest tenth and determining the mode of each row in a frame. These 
% modes can be used to identify the displacements that occur the most which
% which we will assume to be the case for displacement caused by patterning
% error, and not by a cell seeded on top of the gel.
%
% We then take the average value of x and y displacement of every row 
% member yielding the x AND y mode values and substract those values from 
% the starting x and y positions of the trajectories in the current test 
% frame.
%
% By shifting the start position of the vectors, we can remove the error in
% displacement caused by patterning irregularities

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.1 Identifying Deviations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if attemptRowAdjust == 1
numTrajInRows = max(rows(:));
numRowElements = size(rows,2);
bookX = zeros(numRowElements+100,max(book1(9,1,:)),totalNumFrames);
bookY = zeros(numRowElements+1,max(book1(9,1,:)),totalNumFrames);

    for f=1:totalNumFrames
       for i=1:numTrajInRows
            j=book1(9,f,i);
            if max(bookX(:,j,f))==0
                m=1;
            end
            k=book1(7,f,i);
            l=book1(8,f,i);
            numRowElements = size(rows,2);
            bookX((numRowElements+1-m),j,f)= k;
            bookY((numRowElements+1-m),j,f)= l;
            m=m+1;
       end
    end


    %filter bookX and bookY to remove noise caused by particle tracker and find
    %the mode error from patterning
    fBookX = bookX .* (abs(bookX)>trackErrorThreshold);
    fBookY = bookY .* (abs(bookY)>trackErrorThreshold);
    numRows = max(nonzeros(book1(9,:,:)));
    for f = 1:totalNumFrames  
       for r = 1:numRows
           if max(abs(fBookX(:,r,f)))>0
           fBookX(1,r,f) = mode(nonzeros(fBookX(:,r,f)));
           else
           fBookX(1,r,f) = 0;
           end
       end
    end

    for f = 1:totalNumFrames  
       for r = 1:numRows
           if max(abs(fBookY(:,r,f)))>0
           fBookY(1,r,f) = mode(nonzeros(fBookY(:,r,f)));
           else
           fBookY(1,r,f) = 0;
           end
       end
    end


    %Transcribe modes to book1 from the fBooks
    for j = 1:numTrajInRows
        for f = 1:totalNumFrames
            i = book1(9,f,j);
            book1(13,f,j) = fBookX(1,i,f);
            book1(14,f,j) = fBookY(1,i,f);
        end
    end
end
%Create new start coordinates
for i = 1:numTraj
    for f = 1:totalNumFrames
        if abs(book1(5,f,i)) < trackErrorThreshold
            book1(15,f,i) = book1(3,f,i);
            book1(17,f,i) = book1(5,f,i);
            book1(21,f,i) = book1(1,f,i);
        else
            book1(15,f,i) = book1(3,f,i);%+book1(13,f,i);
            book1(17,f,i) = book1(5,f,i);%-book1(13,f,i);
            book1(21,f,i) = book1(1,f,i);%-book1(13,f,i);
        end
    end
        if i == 1 || i == 5000 || i == 10000 || i == 15000
        progress = 'Section 4.1 x'
    end
end
for i = 1:numTraj
    for f = 1:totalNumFrames
        if abs(book1(6,f,i)) < trackErrorThreshold
            book1(16,f,i) = book1(4,f,i);
            book1(18,f,i) = book1(6,f,i);
            book1(22,f,i) = book1(2,f,i);
        else
            book1(16,f,i) = book1(4,f,i);% +book1(14,f,i); add these to track constant y errors
            book1(18,f,i) = book1(6,f,i);% -book1(14,f,i); right now they are not necessary
            book1(22,f,i) = book1(2,f,i);% -book1(14,f,i); 
        end
    end
        if i == 1 || i == 5000 || i == 10000 || i == 15000
        progress = 'Section 4.1 y'
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.2 Storing the Maximum Displacement Values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: Some trajectories do not persist all the way to the top frame, and
%so they would not be visible in the outputs in section 5 unless their
%maximum value is used for those overlays.

for i = 1:numTraj
   
    book1(23,:,i) = book1(17,book1(12,1,i),i);
    book1(24,:,i) = book1(18,book1(12,1,i),i);
    book1(25,:,i) = book1(21,book1(12,1,i),i);
    book1(26,:,i) = book1(22,book1(12,1,i),i);
    book1(29,:,i) = (book1(23,:,i).^2 + book1(24,:,i).^2).^0.5;
    
    if i == 1 || i == 5000 || i == 10000 || i == 15000
    progress = 'Section 4.2'
    end
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.3 Creating a Color Map for Quiver Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

promptColorMap = 'How many colors on color map? Enter an integer and press enter: ';
cMD = input(promptColorMap) %color map divisions
[map,scheme] = brewermap_view(cMD)
cMD = size(map,1)
divSize = (max(max(book1(19,:,:))))/cMD;
cMDBook = zeros(numIndices,totalNumFrames,numTraj,cMD); %for creating colormaps of deformations in current frame
cMDBook2 = zeros(numIndices,totalNumFrames,numTraj,cMD); %for creating colormaps of largest deformation all in one frame
maskArray = zeros(numIndices,totalNumFrames,numTraj,1);
for i = 1:cMD
     cMDBook(:,:,:,i) = book1(:,:,:);     
     maskArray(1,:,:,1) = book1(19,:,:) < (divSize*i) & book1(19,:,:) > (divSize*(i-1));
     for j = 1:numIndices         
     cMDBook(j,:,:,i) = cMDBook(j,:,:,i).*maskArray(1,:,:,1);
     end
end
progress = 'done cMDBook1'
for i = 1:cMD
     cMDBook2(:,:,:,i) = book1(:,:,:);
     maskArray(1,:,:,1) = book1(29,:,:) < (divSize*i) & book1(29,:,:) > (divSize*(i-1));
     for j = 1:numIndices         
     cMDBook2(j,:,:,i) = cMDBook2(j,:,:,i).*maskArray(1,:,:,1);
     end
end
progress = 'done cMDBook2'
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.4a Using Intensity Values to Extract Z-Information%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: Calculates average distance between PSFs in a single pillar.
%Outputs stored as variable "yDataDiffAverage"
if ismember(10,outputs) == 1
close all

intProfiles = zeros(totalNumFrames,numTraj);

sortedArray = zeros(numTraj,3);
for i = 1:numTraj
    sortedArray(i,1) = i;
    sortedArray(i,2) = book1(3,1,i);
    sortedArray(i,3) = book1(4,1,i);
end

[sortedArray,sortedArrayIndices] = sortrows(sortedArray,[3, 2]);

for i = 1:(numTraj)
    intProfiles(book1(11,1,sortedArray(i,1)):book1(12,1,sortedArray(i,1)),i) = book1(30,book1(11,1,sortedArray(i,1)):book1(12,1,sortedArray(i,1)),sortedArray(i,1));    
end
%image(intProfiles,'cdatamapping','scaled'); 
%colormap(gray);
%imwrite(intProfilesImage,'intProfilesImage.tiff');
close;
intProfilesPS = phasesym(intProfiles,'minWaveLength',3,'mult',1.5);
intProfilesMax = (localmax(intProfilesPS'))' > 0;
intProfilesClosed1 = imclose(intProfilesMax, ones(1,20));
intProfilesCleaned = bwareaopen(intProfilesClosed1,200);
intProfilesClosed2 = imclose(intProfilesCleaned, ones(1,1000));

intProfilesLabeled = bwlabel(intProfilesClosed2);

intProfilesFitCoefs = zeros(5,5);
yData = zeros(numTraj,5);
yDataDiffAverage = zeros(1,5);
intProfilesFits = figure;
image(intProfiles,'cdatamapping','scaled'); colormap(gray);
hold on
for i = 1:max(max(intProfilesLabeled))
    [intProfilesCurrentLabely,intProfilesCurrentLabelx] = find((intProfilesLabeled == i).*intProfilesCleaned);
    intProfilesFitCoefs(1:2,i) = polyfit(intProfilesCurrentLabelx, intProfilesCurrentLabely, 1);
    xData = linspace(1,numTraj,numTraj);
    yData(1:numTraj,i) = intProfilesFitCoefs(1,i)*(xData) + intProfilesFitCoefs(2,i);
    plot(xData,yData(1:numTraj,i))
    hold on
    plot(intProfilesCurrentLabelx, intProfilesCurrentLabely,'*','MarkerSize',5);
    hold on
end
savefile = [trajPath '\Intensity Profile Linear Fits.tif'];
export_fig(intProfilesFits,savefile);
hold off

intProfilesRaw = figure;
image(intProfiles,'cdatamapping','scaled'); colormap(gray);
savefile = [trajPath '\Intensity Profiles 1 Raw.tif'];
export_fig(intProfilesRaw,savefile);

intProfilesPSfig = figure;
image(intProfilesPS,'cdatamapping','scaled'); colormap(gray);
savefile = [trajPath '\Intensity Profiles 2 Phase Sym.tif'];
export_fig(intProfilesPSfig,savefile);

intProfilesMaxfig = figure;
image(intProfilesMax,'cdatamapping','scaled'); colormap(gray);
savefile = [trajPath '\Intensity Profiles 3 Max.tif'];
export_fig(intProfilesMaxfig,savefile);

intProfilesClosedfig1 = figure;
image(intProfilesClosed1,'cdatamapping','scaled'); colormap(gray);
savefile = [trajPath '\Intensity Profiles 4 Closed.tif'];
export_fig(intProfilesClosedfig1,savefile);

intProfilesCleanedfig = figure;
image(intProfilesCleaned,'cdatamapping','scaled'); colormap(gray);
savefile = [trajPath '\Intensity Profiles 5 Cleaned.tif'];
export_fig(intProfilesCleanedfig,savefile);

intProfilesClosedfig2 = figure;
image(intProfilesClosed2,'cdatamapping','scaled'); colormap(gray);
savefile = [trajPath '\Intensity Profiles 6 Closed.tif'];
export_fig(intProfilesClosedfig2,savefile);

intProfilesLabeledfig = figure;
image(intProfilesLabeled,'cdatamapping','scaled'); colormap(gray);
savefile = [trajPath '\Intensity Profiles 7 Labeled.tif'];
export_fig(intProfilesLabeledfig,savefile);

for i = 2:max(max(intProfilesLabeled))
    yDataDiffAverage(1,i) = mean(yData(:,i-1)-yData(:,i));
end





% close;
% intProfilesPSM = phasesymmono(intProfiles,'minWaveLength',2,'mult',1.5);
% intProfilesPS = phasesym(intProfiles,'minWaveLength',2,'mult',1.5);
% intProfilesPCM = phasecongmono(intProfilesPSM);
% [intProfilesPC3,trash, orient] = phasecong3(intProfilesPS);
% [intProfilesNMS, intProfilesNMSPos] = nonmaxsup(intProfilesPS,intProfilesOrient,1.3);
% intProfilesNMS2 = medfilt2(intProfilesNMS, [3 3]) > 0;
% intProfilesNMS3 = imclose(intProfilesNMS2, ones(1,200));
% intProfilesNMSM = nonmaxsuppts(intProfilesPCM);
% 
% image(intProfilesNMS3,'cdatamapping','scaled'); colormap(gray);

%imwrite(intProfilesImage,'intProfilesImage.tiff');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.4b Using Intensity Values to Extract Z-Information%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(9,outputs) == 1
cMDBook2(cMDBook2 == 0) = NaN;
% plottedProfiles = figure
% for i = 1:numTraj
% plot(linspace(1,54,54),cMDBook2(30,:,i,18),'.','MarkerSize',5)
% hold on
% end

trajOfInterest = 1;
winSize = 50
if trajOfInterest > winSize && trajOfInterest < (numTraj-winSize) 
    lowerLim = trajOfInterest - winSize;
    upperLim = trajOfInterest + winSize;
elseif trajOfInterest < winSize
    lowerLim = trajOfInterest-(trajOfInterest-1);
    upperLim = trajOfInterest*2 + winSize;
elseif trajOfInterest > (numTraj-winSize)
    lowerLim = trajOfInterest - winSize - (numTraj-trajOfInterest);
    upperLim = trajOfInterest + (trajOfInterest-1);
end

plottedProfiles = figure
plot(linspace(1,totalNumFrames,totalNumFrames),book1(30,:,trajOfInterest),'.','MarkerSize',5,'Color',[0 0 0])
hold on
findpeaks((book1(30,:,trajOfInterest)))

currentAverages = zeros(totalNumFrames,1);
 
for i = 1:totalNumFrames
errorbar(i,mean(cMDBook2(30,i,lowerLim:upperLim,1),'omitnan'),std(cMDBook2(30,i,:,1),'omitnan'),'.','MarkerSize',5,'Color','r')
hold on
currentAverages(i,1) = mean(cMDBook2(30,i,lowerLim:upperLim,1),'omitnan');
end
currentAverages(isnan(currentAverages))=0;
findpeaks(currentAverages(:,1))

hold on
xVals = linspace(1,totalNumFrames,totalNumFrames);
fitWeights = zeros(totalNumFrames,1);
fitWeights(:,1) = xVals(1,:)+1000;
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[4000,25000,0,0,0,0,0],...
               'Upper',[6000,35000,15,.1,1,3.14,15],...
               'StartPoint',[5500,29000,11.5,.015,0.1731,2.2,6.712]);
f = fittype('pillarFit2(x,S,iS,iP,P,Pg,Mo,Mf)','options',s);
[fitCoef,fitQual] = fit(xVals',currentAverages,f,'weight',fitWeights); 
fitCoef2 = coeffvalues(fitCoef);
yVals = pillarFit(xVals,fitCoef2(1,5),fitCoef2(1,7),fitCoef2(1,6),fitCoef2(1,3),fitCoef2(1,4),fitCoef2(1,2),fitCoef2(1,1));
plot(xVals,yVals,'MarkerSize',20,'Color','bl')
findpeaks(yVals(1,:))


% folderName = strcat( 'Intensity Average Profiles with ',num2str(cMD),' Divisions');
% mkdir(blackPath,folderName)
% for j = 1:cMD
% plottedProfilesErrorBars = figure
% for i = 1:totalNumFrames
% errorbar(i,mean(cMDBook2(30,i,:,j),'omitnan'),std(cMDBook2(30,i,:,j),'omitnan'),'.','MarkerSize',5)
% hold on
% end
% 
% axis([0 60 0 40000])
% savefile = [blackPath '\' folderName '\Displacement Magnitude Division ' num2str(j) '.tif'];
% export_fig(plottedProfilesErrorBars,savefile);
% 
% close
% end
end
%%


%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%5.)IMAGE OUTPUTS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%5.1 Drawing Zero-State Displacement Fields%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismember(1,outputs) == 1
for f = 1:totalNumFrames        % number of z-slices
    res= {152,154};
    data = figure('units','pixels','outerposition',[0 0 resX resY]);
    h = quiver(book1(15,f,:),resY-book1(16,f,:),book1(17,f,:),-book1(18,f,:),0);
    hU = get(h,'UData');
    hV = get(h,'VData');
    set(h,'UData',scaleOutputVectors*hU,'VData',scaleOutputVectors*hV)
    axis(gca,'equal','tight')
    savefile = sprintf('Trajectory between 1st Frame and Frame %u.tif',f);
    print(savefile,'-dtiff','-r300')
    close
end
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%5.2 Drawing Zero-State Displacement Fields on Black Background%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismember(2,outputs) == 1
    folderName = strcat('Black Image_',scheme,' Quiver Overlays');
    mkdir(blackPath,folderName)
    for f = 1:totalNumFrames        % number of z-slices
        blackOverlay = figure('Position',[0 0 1000 1000]);
        imshow(e,[])
        hold on
%       quiver(book1(15,f,:),book1(16,f,:),book1(17,f,:),book1(18,f,:),0,'g');
        for i = 1:cMD
        quiver(cMDBook(15,f,:,i),cMDBook(16,f,:,i),cMDBook(17,f,:,i),cMDBook(18,f,:,i),0,'color',[map(i,1:3)]);
        hold on
        end
        hold off
        savefile = [blackPath '\' folderName '\Black Background Overlay' scheme ' ' num2str(f) '.tif'];
        export_fig(blackOverlay,savefile,'-native');
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
folderName = strcat( 'Centroid_',scheme,' Black Image Overlays');
mkdir(blackPath, folderName)
for f = 1:totalNumFrames
centroidsOnly = figure('Position',[0 0 1000 1000]);
imshow(e,[])
hold on
        for i = 1:cMD
        
        xTemp = squeeze(squeeze(cMDBook(1,f,:,i)));
        xTemp = xTemp';
        yTemp = squeeze(squeeze(cMDBook(2,f,:,i)));
        yTemp = yTemp';
        plot(xTemp,yTemp,'.','MarkerSize',3,'Color',[map(i,1:3)]);
        %quiver(cMDBook(15,f,:,i),cMDBook(16,f,:,i),cMDBook(17,f,:,i),cMDBook(18,f,:,i),0,'color',[red green blue]);
        hold on
        end
% xTemp = squeeze(book1(21,f,:));
% yTemp = squeeze(book1(22,f,:));
% plot(xTemp,yTemp,'w.','MarkerSize',3);
hold on
hold off
savefile = [blackPath '\' folderName '\Centroids on Frame ' scheme ' ' num2str(f) '.tif'];
export_fig(centroidsOnly,savefile,'-native');
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
imshow(d,[])
hold on
xTemp = squeeze(book1(21,1,:));
yTemp = squeeze(book1(22,1,:));
plot(xTemp,yTemp,'r.','MarkerSize',10);
hold on
if ismember(7,outputs) == 0
for f = 1:totalNumFrames
    xTemp = squeeze(book1(21,f,:));
    yTemp = squeeze(book1(22,f,:));
    plot(xTemp,yTemp,'b.','MarkerSize',3);
    hold on
end
end
hold on
if ismember(8,outputs) == 0
    for i = 1:cMD
        quiver(cMDBook2(27,totalNumFrames,:,i),cMDBook2(28,totalNumFrames,:,i),cMDBook2(23,totalNumFrames,:,i),cMDBook2(24,totalNumFrames,:,i),0,'color',[map(i,1:3)]);
        hold on
    end
%quiver(book1(3,1,:),book1(4,1,:),book1(23,totalNumFrames,:),book1(24,totalNumFrames,:),0,'g');
end
hold off
savefile = [blackPath '\Debug Image ' scheme '.tif'];
export_fig(debugImage,savefile,'-native');
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%5.5 Plotting Transmitted Quiver Overlays v1%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: 5.5 Saves files, 5.6 does not.

if ismember(5,outputs) == 1
trajOverlay = figure('Position',[0 0 1000 1000]);
imshow(c,[])
hold on
 for i = 1:cMD

        quiver(cMDBook2(27,totalNumFrames,:,i),cMDBook2(28,totalNumFrames,:,i),cMDBook2(23,totalNumFrames,:,i),cMDBook2(24,totalNumFrames,:,i),(i^1.3)/((cMD/2)^1.3),'color',[map(i,1:3)]);%(i^2)/((cMD/1.5)^2)
        hold on
    end
%quiver(book1(27,1,:),book1(28,1,:),book1(23,totalNumFrames,:),book1(24,totalNumFrames,:),0,'g');
hold off
savefile = [trajPath '\Fluorescent Overlay ' scheme '.tif'];
export_fig(trajOverlay,savefile,'-native');
end

%%
if ismember(5,outputs) == 1
transmittedOverlay = figure('Position',[0 0 1000 1000]);
imshow(d,[])
hold on
 for i = 1:cMD
        quiver(cMDBook2(27,totalNumFrames,:,i),cMDBook2(28,totalNumFrames,:,i),cMDBook2(23,totalNumFrames,:,i),cMDBook2(24,totalNumFrames,:,i),(i^1.3)/((cMD/2)^1.3),'color',[map(i,1:3)]);
        hold on
    end
%quiver(book1(27,1,:),book1(28,1,:),book1(23,totalNumFrames,:),book1(24,totalNumFrames,:),0,'g');
hold off
savefile = [cellPath '\Transmitted Overlay ' scheme '.tif'];
export_fig(transmittedOverlay,savefile,'-native');
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%5.6 Plot Overlay on Trajectories and Transmitted Images v2%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: 5.5 Saves files, 5.6 does not.

if ismember(6,outputs) == 1
figure('Position',[0 0 1000 1000])
imshow(c,[])
hold on
 for i = 1:cMD
        if (i/cMD) >= .5
            blue = 0;
            green = 1-(((i/cMD)-0.5)/(.5));
            red = 1 - green ;
        elseif (i/cMD) < .5
            blue = ((0.5-(i/cMD))/0.5);
            green = 1 - blue ;
            red = 0;
        end
        quiver(cMDBook2(27,totalNumFrames,:,i),cMDBook2(28,totalNumFrames,:,i),cMDBook2(23,totalNumFrames,:,i),cMDBook2(24,totalNumFrames,:,i),0,'color',[map(i,1:3)]);
        hold on
    end
%quiver(book1(27,1,:),book1(28,1,:),book1(23,totalNumFrames,:),book1(24,totalNumFrames,:),0,'g');
hold off

figure('Position',[0 0 1000 1000])
imshow(d,[])
hold on
 for i = 1:cMD
        if (i/cMD) >= .5
            blue = 0;
            green = 1-(((i/cMD)-0.5)/(.5));
            red = 1 - green ;
        elseif (i/cMD) < .5
            blue = ((0.5-(i/cMD))/0.5);
            green = 1 - blue ;
            red = 0;
        end
        quiver(cMDBook2(27,totalNumFrames,:,i),cMDBook2(28,totalNumFrames,:,i),cMDBook2(23,totalNumFrames,:,i),cMDBook2(24,totalNumFrames,:,i),0,'color',[map(i,1:3)]);
        hold on
    end
%quiver(book1(27,1,:),book1(28,1,:),book1(23,totalNumFrames,:),book1(24,totalNumFrames,:),0,'g');
hold off
end

%%

