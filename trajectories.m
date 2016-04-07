%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converting Particle Tracking Data to Displacement Data from Zero State  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
%1.) Sorting the raw data from excel file
%--------------------------------------------------------------------------
close all; clear; clc;

% for u(ser)
[name,path] = uigetfile('*.xlsx');
file = [path,name];

% % for me
% file = 'F:\trajectories.xlsx';
[num,txt,raw] = xlsread(file);
resX = 152;
resY = 274;
scaleOutputVectors = 1;
trackErrorThreshold = .3;
patternErrorThreshold = .5;
numIndices = 20;
traj = num(:,2);
numTraj = max(traj); % number of objects
totalNumFrames = max(num(:,3)) + 1; % Maximum number of frames observable for any one object
book1 = zeros(numIndices,totalNumFrames,numTraj);
for i = 1:numTraj
    %Here we build a book of pages (3D array) with the data for a single
    %object/trajectory per page. Most of the code is to ensure that each
    %page is the same size matrix as the next.
    
    tempObj = num(traj==i,:); %stores data to 'tempObj' pertaining to all frames of one object/trajectory
    numFrames = size(tempObj,1); % number of frames that the current object appears in.
    startFrame = min(tempObj(:,3)); % the first frame that an object appears in (tracking software starts at 0)
    endFrame = (startFrame+numFrames); % the last frame that an object appears in
    book1(11,:,i) = startFrame+1;
    book1(12,:,i) = endFrame;
    if startFrame > 0  % this statement fills in the upper portion of obj matrix with zeros if the first frame is not 0  
        obj = [zeros(startFrame,12);tempObj];
    else
        if numFrames < totalNumFrames
        obj = [tempObj;zeros(totalNumFrames-numFrames,12)];
        else
        if startFrame == 0
        obj = tempObj; %this statement applys to objects that start in frame 0
        end
        end
    end
    if i == 1
        book2 = zeros(totalNumFrames, 12); %'book1' is our book (3D array) of object data
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
%2.) Building book1 from the raw data
%--------------------------------------------------------------------------
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
%3.) Associating objects to rows by identifying nearest right neighbor 
%--------------------------------------------------------------------------
objAll = [x_0(:,1),resY-y_0(:,1)];
rows = zeros(2,1);
row = 1;
col = 1;
% distances = cell(numObj,3);
%for i = 1:21
for i = 1:numTraj
    if book1(11,1,i) ~= 1
        i=i+1;
    else    
    thisObj = [x_0(i,1),resY-y_0(i,1)];
    objCount = (1:numTraj)';
    distances = sqrt(sum(bsxfun(@minus,objAll,thisObj).^2,2));
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
    if isempty (rows)
        rows(row,col) = i;
        book1(10,:,i)= rightNeighbor(3);
        book1(9,:,i) = row;
    else
        if isempty(rightNeighbor)   %if a right neighbor doesn't exist for the current trajectory (i.e. if it is on the edge of the frame)
            col = col + 1;           %adjust the current column to the next slot(to the right), and
            rows(row,col) = i;      %place the current trajectory in that slot
            book1(9,:,i) = row;
            row = row + 1;          %then start a new row for the next trajectory in the list (because we have reached an edge)
            col = 1;                %start at the beginning of the new row for the next trajectory
            book1(10,:,i)= 0;
            
        else                                                     %if we are NOT at the end of the row, we assume first that we MAY NOT be in the correct row, 
                                                                 %we check the current trajectory's right neighbor's row to see if it exists 
            rowCheck = ismember(rightNeighbor(3),rows);          %is the right neighbor to current object already assigned a row? store as 'rowCheck'
            if rowCheck == 1                                     %if yes then
                [trow,tcol] = find(rows == rightNeighbor(3));    %what row does the right neighbor belong to? Store as 'trow'
                emptyColumn = size(rows,2)+1;                    %Add one to the width of the 'rows' array to have a slot which will necessarily be empty, stored as 'emptyColumn'
                if rows(trow,emptyColumn-1) == 0                 %If there is an empty value to the left of 'emptyColumn' in 'trow' then,
                rows(trow,emptyColumn-1) = i;                    %store trajectory number in that empty slot
                else                                             %otherwise,
                rows(trow,emptyColumn) = i;                      %place the trajectory number into a new empty column
                end
                book1(10,:,i)= rightNeighbor(3);
                book1(9,:,i) = trow;
            else
                col = col + 1;                                   %if the right neighbor is not on the list of rows yet, then move one column right and
                rows(row,col) = i;                               %add the current trajectory to the empty slot
                book1(10,:,i)= rightNeighbor(3);                 %add rightNeighbor to index book1 (index#10)
                book1(9,:,i) = row;
            end
        end
    end
    
    end
end

%%
%--------------------------------------------------------------------------
%4.)Identifying deviations from zero state in later frames
%--------------------------------------------------------------------------
%Now that we have a list of rows of trajectories we need to find the
%error displacement for members of a particular row caused by irregularities in the patterning process. 
%
%We are measuring from 
%the zero frame to
%the current frame. 
%
%We can do this by rounding the x and y displacement of each value to the nearest
%tenth and determining the mode of each row in a frame. These modes can be
%used to identify the displacements that occur the most which which we will assume
%to be the case for displacement caused by patterning error, and not by a
%cell seeded on top of the gel.
%
%We then take the average value of x and y displacement of every row member
%yielding the x AND y mode values and substract those values from the starting x and
%y positions of the trajectories in the current test frame.
%
%By shifting the start position of the vectors, we can remove the error in
%displacement caused by patterning irregularities

numTrajInRows = max(rows(:));
rowElements = find(book1(9,1,:)==1);
numRowElements = numel(rowElements);
bookX = zeros(numRowElements+1,max(book1(9,1,:)),totalNumFrames);
bookY = zeros(numRowElements+1,max(book1(9,1,:)),totalNumFrames);
for f=1:totalNumFrames
    m=2;
    for i=1:numTrajInRows
    j=book1(9,f,i);
    k=book1(7,f,i);
    l=book1(8,f,i);
    rowElements = find(book1(9,1,:)==j);
    numRowElements = numel(rowElements);
    bookX(m-((j-1)*numRowElements),j,f)= k;
    bookY(m-((j-1)*numRowElements),j,f)= l;
    m=m+1;
    end
end


%filter bookX and bookY to remove noise caused by particle tracker and find
%the mode error from patterning
fBookX = bookX .* (abs(bookX)>trackErrorThreshold);
fBookY = bookY .* (abs(bookX)>trackErrorThreshold);
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
        else
            book1(15,f,i) = book1(3,f,i)-book1(13,f,i);
            book1(17,f,i) = book1(5,f,i)-book1(13,f,i);
        end
    end
end
for i = 1:numTraj
    for f = 1:totalNumFrames
        if abs(book1(6,f,i)) < trackErrorThreshold
            book1(16,f,i) = book1(4,f,i);
            book1(18,f,i) = book1(6,f,i);
        else
            book1(16,f,i) = book1(4,f,i)-book1(14,f,i);
            book1(18,f,i) = book1(6,f,i)-book1(14,f,i);
        end
    end
end

%%
%--------------------------------------------------------------------------
%Drawing zero state displacement fields
%--------------------------------------------------------------------------
% for f = 1:totalNumFrames        % number of z-slices
%     res= {152,154};
%     data = figure('units','pixels','outerposition',[0 0 resX resY]);
%     h = quiver(book1(15,f,:),resY-book1(16,f,:),book1(17,f,:),-book1(18,f,:),0);
%     hU = get(h,'UData');
%     hV = get(h,'VData');
%     set(h,'UData',scaleOutputVectors*hU,'VData',scaleOutputVectors*hV)
%     axis(gca,'tight')
%     savefile = sprintf('Trajectory between 1st Frame and Frame %u.tif',f);
%     print(savefile,'-dtiff','-r300')
%     close
% end

%%
[a,b] = uigetfile('*.tif','Select trajectories image');
c = imread([b,a]);
figure
imshow(c,[])
hold on
quiver(book1(15,totalNumFrames,:),book1(16,totalNumFrames,:),book1(17,totalNumFrames,:),book1(18,totalNumFrames,:),0,'g');
hold off

[a,b] = uigetfile('*.tif','Select cell overlay image');
d = imread([b,a]);
figure
imshow(d,[])
hold on
quiver(book1(15,totalNumFrames,:),book1(16,totalNumFrames,:),book1(17,totalNumFrames,:),book1(18,totalNumFrames,:),0,'g');
hold off