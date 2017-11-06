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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.1 Identifying Deviations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            book1(3,f,i) = book1(3,f,i);
            book1(5,f,i) = book1(5,f,i);
            book1(1,f,i) = book1(1,f,i);
        else
            book1(3,f,i) = book1(3,f,i);%+book1(13,f,i);
            book1(5,f,i) = book1(5,f,i);%-book1(13,f,i);
            book1(1,f,i) = book1(1,f,i);%-book1(13,f,i);
        end
    end
        if i == 1 || i == 5000 || i == 10000 || i == 15000
        progress = 'Section 4.1 x'
    end
end
for i = 1:numTraj
    for f = 1:totalNumFrames
        if abs(book1(6,f,i)) < trackErrorThreshold
            book1(4,f,i) = book1(4,f,i);
            book1(6,f,i) = book1(6,f,i);
            book1(2,f,i) = book1(2,f,i);
        else
            book1(4,f,i) = book1(4,f,i);% +book1(14,f,i); add these to track constant y errors
            book1(6,f,i) = book1(6,f,i);% -book1(14,f,i); right now they are not necessary
            book1(2,f,i) = book1(2,f,i);% -book1(14,f,i); 
        end
    end
        if i == 1 || i == 5000 || i == 10000 || i == 15000
        progress = 'Section 4.1 y'
    end
end
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
%%%%%%%%%5.6 Plot Overlay on Trajectories and Transmitted Images v2%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: 5.5 Saves files, 5.6 does not.

if ismember(6,outputs) == 1
    figure('Position',[0 0 1000 1000])
    imshow(imageFluor,[])
    hold on
    for i = 1:cmD
        if (i/cmD) >= .5
            blue = 0;
            green = 1-(((i/cmD)-0.5)/(.5));
            red = 1 - green ;
        elseif (i/cmD) < .5
            blue = ((0.5-(i/cmD))/0.5);
            green = 1 - blue ;
            red = 0;
        end
        quiver(cm2(27,totalNumFrames,:,i),cm2(28,totalNumFrames,:,i),cm2(23,totalNumFrames,:,i),cm2(24,totalNumFrames,:,i),0,'color',[colorMap(i,1:3)]);
        hold on
    end
    %quiver(book1(10,1,:),book1(11,1,:),book1(12,totalNumFrames,:),book1(13,totalNumFrames,:),0,'g');
    hold off
    
    figure('Position',[0 0 1000 1000])
    imshow(imageTrans,[])
    hold on
    for i = 1:cmD
        if (i/cmD) >= .5
            blue = 0;
            green = 1-(((i/cmD)-0.5)/(.5));
            red = 1 - green ;
        elseif (i/cmD) < .5
            blue = ((0.5-(i/cmD))/0.5);
            green = 1 - blue ;
            red = 0;
        end
        quiver(cm2(27,totalNumFrames,:,i),cm2(28,totalNumFrames,:,i),cm2(23,totalNumFrames,:,i),cm2(24,totalNumFrames,:,i),0,'color',[colorMap(i,1:3)]);
        hold on
    end
    %quiver(book1(10,1,:),book1(11,1,:),book1(12,totalNumFrames,:),book1(13,totalNumFrames,:),0,'g');
    hold off
end