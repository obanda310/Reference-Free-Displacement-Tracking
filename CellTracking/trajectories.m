close all; clear; clc;
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

[trajFile,trajPath] = uigetfile('*.tif','Select trajectories image');
c = imread([trajPath,trajFile]);

[cellFile,cellPath] = uigetfile('*.tif','Select cell overlay image');
d = imread([cellPath,cellFile]);

[blackFile,blackPath] = uigetfile('*.tif','Select black image of correct size');
e = imread([blackPath,blackFile]);

outputs = OutputSelector()

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.2 Input Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[resY,resX] = size(c); %for Section 5
scaleOutputVectors = 1; %for Section 5.1
trackErrorThreshold = .2; %for Section 4
numIndices = 26; %for Section 1 and 2 (number of elements in book1)
numTraj = max(num(:,2)); %Number of Trajectories
totalNumFrames = max(num(:,3)) + 1; %Maximum number of frames observable for any one object
book1 = zeros(numIndices,totalNumFrames,numTraj); %Creates book1

%%

for i = 1:numTraj
    %Here we build a book of pages (3D array) with the data for a single
    %object/trajectory per page. Most of the code is to ensure that each
    %page is the same size matrix as the next.
    
    % stores data to 'tempObj' pertaining to all frames of one
    % object/trajectory
    tempObj = num(num(:,2)==i,:);
    % number of frames that the current object appears in.
    numFrames = size(tempObj,1);
    % the first frame that an object appears in (tracking software starts 
    % at 0
    startFrame = min(tempObj(:,3));
    % the last frame that an object appears in
    endFrame = (startFrame+numFrames);
    book1(11,:,i) = startFrame+1;
    book1(12,:,i) = endFrame;
    % this fills in the upper portion of obj matrix with zeros if the first
    % frame is not 0
    if startFrame > 0 && endFrame < totalNumFrames
        obj = [zeros(startFrame,12);tempObj;zeros(totalNumFrames-endFrame,12)];
    else
    if endFrame < totalNumFrames
        obj = [tempObj;zeros(totalNumFrames-endFrame,12)];
    else
    if startFrame > 0  
        obj = [zeros(startFrame,12);tempObj];
    else
        if numFrames < totalNumFrames
        obj = [tempObj;zeros(totalNumFrames-numFrames,12)];
        else
        if startFrame == 0
        obj = tempObj; %this applies to objects that start in frame 0
        end
        end
    end
    end
    end
    if i == 1
        %'book1' is our book (3D array) of object data
        book2 = zeros(totalNumFrames, 12);
        book2(:,:) = obj;
    else
        book2 = cat(3,book2,obj); %adding pages of 'obj' data to the book.
    end
   
 %Here we use the temporary data stored in 'obj' to build 2 matrices to
 %describe the x and y positions in frame 1 for each object, as well as the
 %displacement from those matrices occurring in successive frames.
    for j = startFrame+1:endFrame
        x_0(i,j) = obj(startFrame+1,4);
        y_0(i,j) = obj(startFrame+1,5);
        x_diff(i,j) = obj(j,4) - obj(startFrame+1,4);
        y_diff(i,j) = obj(j,5) - obj(startFrame+1,5);
    end
end
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
%19 = Currently empty
%20 = Nearest Left neighbor Traj to current Traj
%21 = Mode Corrected X Coordinate of Traj in Current Frame
%22 = Mode Corrected Y Coordinate of Traj in Current Frame
%23 = Value of book1 index 17 (see above) in Last Frame that Traj Appears
%24 = Value of book1 index 18 (see above) in Last Frame that Traj Appears
%25 = Value of book1 index 21 (see above) in Last Frame that Traj Appears
%26 = Value of book1 index 22 (see above) in Last Frame that Traj Appears

  for i = 1:numTraj
        for j = 1:totalNumFrames
            book1(1,j,i) = book2(j,4,i);
            book1(2,j,i) = book2(j,5,i);
            book1(3,j,i) = x_0(i,j);
            book1(4,j,i) = y_0(i,j);
            book1(5,j,i) = x_diff(i,j);
            book1(6,j,i) = y_diff(i,j);
            book1(7,j,i) = round(book1(5,j,i),1);
            book1(8,j,i) = round(book1(6,j,i),1);
        end
  end

%%


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 3.) Assigning Objects to Rows by Based on Neighbor Particles
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%3.1 Identifying Nearest Left and Right Neighbors%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

objAll = [x_0(:,1),resY-y_0(:,1)];
rows = zeros(2,1);
row = 1;
col = 1;
for i = 1:numTraj
    if book1(11,1,i) ~= 1
        i=i+1;
    else    
        thisObj = [x_0(i,1),resY-y_0(i,1)];
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
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.2 Storing the Maximum Displacement Values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: Some trajectories do not persist all the way to the top frame, and
%so they would not be visible in the outputs in section 5 unless their
%maximum value is used for those overlays.

for i = 1:numTraj
    ind17 = find(book1(17,:,i));
    ind17 = ind17(end);
    book1(23,:,i) = book1(17,ind17,i);
    ind18 = find(book1(18,:,i));
    ind18 = ind18(end);
    book1(24,:,i) = book1(18,ind18,i);
    ind21 = find(book1(21,:,i));
    ind21 = ind21(end);
    book1(25,:,i) = book1(21,ind21,i);
    ind22 = find(book1(22,:,i));
    ind22 = ind22(end);
    book1(26,:,i) = book1(22,ind22,i);
end

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
    mkdir(blackPath,'Black Image Quiver Overlays')
    for f = 1:totalNumFrames        % number of z-slices
        blackOverlay = figure
        imshow(e,[])
        hold on
        quiver(book1(15,f,:),book1(16,f,:),book1(17,f,:),book1(18,f,:),0,'g');
        hold off
        savefile = [blackPath '\Black Image Quiver Overlays\Black Background Overlay' num2str(f) '.tif'];
        export_fig(blackOverlay,savefile);
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
mkdir(blackPath,'Centroid Black Image Overlays')
for f = 1:totalNumFrames
centroidsOnly = figure;
imshow(e,[])
hold on
xTemp = squeeze(book1(21,f,:));
yTemp = squeeze(book1(22,f,:));
plot(xTemp,yTemp,'w.','MarkerSize',3);
hold on
hold off
savefile = [blackPath '\Centroid Black Image Overlays\Centroids on Frame ' num2str(f) '.tif'];
export_fig(centroidsOnly,savefile);
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
debugImage = figure('units','pixels','outerposition',[0 0 resX resY]);
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
quiver(book1(3,1,:),book1(4,1,:),book1(23,totalNumFrames,:),book1(24,totalNumFrames,:),0,'g');
end
hold off
savefile = [blackPath '\Debug Image.tif'];
export_fig(debugImage,savefile);
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%5.5 Plotting Transmitted Quiver Overlays v1%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: 5.5 Saves files, 5.6 does not.

if ismember(5,outputs) == 1
trajOverlay = figure;
imshow(c,[])
hold on
quiver(book1(3,totalNumFrames,:),book1(4,totalNumFrames,:),book1(23,totalNumFrames,:),book1(24,totalNumFrames,:),0,'g');
hold off
savefile = [trajPath '\Trajectories Overlay.tif'];
export_fig(trajOverlay,savefile);

transmittedOverlay = figure;
imshow(d,[])
hold on
quiver(book1(3,1,:),book1(4,1,:),book1(23,totalNumFrames,:),book1(24,totalNumFrames,:),0,'g');
hold off
savefile = [cellPath '\Transmitted Overlay.tif'];
export_fig(transmittedOverlay,savefile);
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%5.6 Plot Overlay on Trajectories and Transmitted Images v2%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: 5.5 Saves files, 5.6 does not.

if ismember(6,outputs) == 1
figure
imshow(c,[])
hold on
quiver(book1(3,1,:),book1(4,1,:),...
    book1(23,totalNumFrames,:),book1(24,totalNumFrames,:),0,'g');
hold off

figure
imshow(d,[])
hold on
quiver(book1(3,1,:),book1(4,1,:),...
    book1(23,totalNumFrames,:),book1(24,totalNumFrames,:),0,'g');
hold off
end

%%
Test(1:30,2) = book1(10,9,1:30);
Test(1:30,1) = book1(20,9,1:30);

