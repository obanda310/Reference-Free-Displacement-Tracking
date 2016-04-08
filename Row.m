classdef Row < Experiment
    properties
        rightNeighbor % trajectory to right of this trajectory
        
        row      % row that trajectory is a member of
        modeRowX   % mode of all round_dx0 trajectory values in this row > threshold
        modeRowY   % mode of all round_dx0 trajectory values in this row > threshold
        newX    % error-adjusted x value
        newY    % error-adjusted y value
        newdx0  % error-adjusted dx0 value
        newdy0  % error-adjusted dx0 value
    end
    methods
        function obj = Row(trajs)
            allIni = {trajs.posIni};
            x0 = cellfun(@(v) v(1), allIni);
            y0 = cellfun(@(v) v(2), allIni);
            allIni = [x0',y0'];
            
            numTrajs = length(trajs);
            count = (1:numTrajs)';
            for i = 1:numTrajs
                if trajs(i).startFrame ~= 1
                    continue
                end
                distances = sqrt(sum(bsxfun(@minus,allIni,trajs(i).posIni).^2,2));
                objData = [count,distances];
                sortDist = sortrows(objData,2);
                % exclude first in sorted distances because that is self-referential
                neighbors = sortDist(2:5,1); % traj ID of neighbors
                nDistance = sortDist(2:5,2); % distance of neighbors
                n = [allIni(neighbors(1,1:4),:),neighbors(1:4)];
                a = bsxfun(@minus,n(:,1:2),trajs(i).posIni);
                b = mean([min(abs(a)),max(abs(a))]);
                goodNeighbors = n(abs(a(:,2)) < b & a(:,1) > 0,:);
                if size(goodNeighbors,1) > 1
                    sortN = sortrows(goodNeighbors,1);
                    rightNeighbor = sortN(1,:);
                else
                    rightNeighbor = goodNeighbors;
                end
            end
        end
    end
end

for i = 1:numTraj
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