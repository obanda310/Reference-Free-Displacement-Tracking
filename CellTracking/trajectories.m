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
[num,txt,raw] = xlsread(file);

[trajFile,trajPath] = uigetfile('*.tif','Select trajectories image');
c = imread([trajPath,trajFile]);
[cellFile,cellPath] = uigetfile('*.tif','Select cell overlay image');
d = imread([cellPath,cellFile]);

[resY,resX] = size(c);

scaleOutputVectors = 1;
trackErrorThreshold = .2;
patternErrorThreshold = .5;
numIndices = 20;
traj = num(:,2);
numTraj = max(traj); % number of objects
% Maximum number of frames observable for any one object
totalNumFrames = max(num(:,3)) + 1;
book1 = zeros(numIndices,totalNumFrames,numTraj);
for i = 1:numTraj
    % Here we build a book of pages (3D array) with the data for a single
    % object/trajectory per page. Most of the code is to ensure that each
    % page is the same size matrix as the next.
    
    % stores data to 'tempObj' pertaining to all frames of one
    % object/trajectory
    tempObj = num(traj==i,:);
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
    goodNeighborsForLeft = n(abs(a(:,2)) < b & a(:,1) < 0,:);
    if size(goodNeighborsForLeft,1) > 1
        sortN = sortrows(goodNeighborsForLeft,1);
        leftNeighbor = sortN(1,:);
    else
        leftNeighbor = goodNeighborsForLeft;
    end
    % if a right neighbor doesn't exist for the current trajectory (i.e. 
    % if it is on the edge of the frame)
        if isempty(rightNeighbor)   
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
            book1(10,:,i)= 0;
        else
            % if we are NOT at the beginning or end of the row, we assume 
            % first that we MAY NOT be in the correct row, we check the 
            % current trajectory's right/left neighbor's row to see if they
            % exist 
            rowCheckRight = ismember(rightNeighbor(3),rows);
            % is the right neighbor to current object already assigned a 
            % row? store as 'rowCheck'
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
                % what row does the right neighbor belong to? Store as trow
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
                book1(10,:,i)= rightNeighbor(3);
                book1(9,:,i) = row;
            end
        end
    end
end
%%
numTrajInRows = max(rows(:));
% for m = 1:size(rows,2)
% for i = 1:numTrajInRows
%     for j = 1:size(rows,1)
%         for k = 1:size(rows,2)
%             if k>1
%             if rows(j,k) == i && rows(j,k-1) == 0
%                 rows(j,k) = 0;
%                 rows(j,k-1) = i;
%                 k = k-2;
%             end
%             end
%         end
%     end
% end
% end
                
                
        
        

%%
%--------------------------------------------------------------------------
%4.)Identifying deviations from zero state in later frames
%--------------------------------------------------------------------------
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


numRowElements = size(rows,2);
bookX = zeros(numRowElements+10,max(book1(9,1,:)),totalNumFrames);
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
            book1(15,f,i) = book1(1,f,i);
            book1(17,f,i) = book1(5,f,i);
        else
            book1(15,f,i) = book1(1,f,i)-book1(13,f,i);
            book1(17,f,i) = book1(5,f,i)-book1(13,f,i);
        end
    end
end
for i = 1:numTraj
    for f = 1:totalNumFrames
        if abs(book1(6,f,i)) < trackErrorThreshold
            book1(16,f,i) = book1(2,f,i);
            book1(18,f,i) = book1(6,f,i);
        else
            book1(16,f,i) = book1(2,f,i);% -book1(14,f,i); add these to track constant y errors
            book1(18,f,i) = book1(6,f,i);% -book1(14,f,i); right now they are not necessary
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
%     axis(gca,'equal','tight')
%    savefile = sprintf('Trajectory between 1st Frame and Frame %u.tif',f);
%     print(savefile,'-dtiff','-r300')
%     close
% end
%%

%--------------------------------------------------------------------------
%Drawing zero state displacement fields on black background
%--------------------------------------------------------------------------
[blackFile,blackPath] = uigetfile('*.tif','Select black image of correct size');
e = imread([blackPath,blackFile]);
for f = 1:totalNumFrames        % number of z-slices
    blackOverlay = figure
    imshow(e,[])
    hold on
    quiver(book1(15,f,:),book1(16,f,:),book1(17,f,:),book1(18,f,:),0,'g');
    hold off
    mkdir(blackPath,'Black Image Overlays')
    savefile = [blackPath '\Black Image Overlays\Black Background Overlay' num2str(f) '.tif'];
    export_fig(blackOverlay,savefile)
    close
end
%%
%--------------------------------------------------------------------------
%Plotting corrected x,y coordinate fields
%--------------------------------------------------------------------------
%data = figure('units','pixels','outerposition',[0 0 resX resY]);
% figure
% imshow(d,[])
% hold on
% for t = 1:numTrajInRows
%     
%     plot(book1(15,1,t),book1(16,1,t),'r.','MarkerSize',10);
%     hold on
%     plot(book1(1,1,t),book1(2,1,t),'g.','MarkerSize',5);
%     hold on
%     for f = 1:totalNumFrames
%     plot(book1(15,f,t),book1(16,f,t),'b.','MarkerSize',3);
%     hold on
% %     plot(book1(1,f,t),book1(2,f,t),'w.','MarkerSize',2);
% %     hold on
%     end
% %     savefile = sprintf('Trajectory between 1st Frame and Frame %u.tif',f);
% %     print(savefile,'-dtiff','-r300')
% %     close
% end
% hold on
% quiver(book1(3,totalNumFrames,:),book1(4,totalNumFrames,:),book1(17,totalNumFrames,:),book1(18,totalNumFrames,:),0,'g');

% %%
% trajOverlay = figure
% imshow(c,[])
% hold on
% quiver(book1(15,totalNumFrames,:),book1(16,totalNumFrames,:),book1(17,totalNumFrames,:),book1(18,totalNumFrames,:),0,'g');
% hold off
% savefile = [trajPath '\Trajectories Overlay.tif'];
% %mkdir('trajPath','Trajectories Overlay')
% export_fig(trajOverlay,savefile)
% 
% transmittedOverlay = figure
% imshow(d,[])
% hold on
% quiver(book1(15,totalNumFrames,:),book1(16,totalNumFrames,:),book1(17,totalNumFrames,:),book1(18,totalNumFrames,:),0,'g');
% hold off
% savefile = [cellPath '\Transmitted Overlay.tif'];
% %mkdir('cellPath','Transmitted Overlay')
% export_fig(transmittedOverlay,savefile)

% figure
% imshow(d,[])
% hold on
% quiver(book1(15,q,:),book1(16,q,:),book1(17,q,:),book1(18,q,:),0,'g');
% hold off
%%=======
figure
imshow(c,[])
hold on
quiver(book1(15,totalNumFrames,:),book1(16,totalNumFrames,:),...
    book1(17,totalNumFrames,:),book1(18,totalNumFrames,:),0,'g');
hold off

figure
imshow(d,[])
hold on
quiver(book1(15,totalNumFrames,:),book1(16,totalNumFrames,:),...
    book1(17,totalNumFrames,:),book1(18,totalNumFrames,:),0,'g');
hold off

