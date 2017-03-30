function [r,rows] = buildRows(r,rowV)
rowV = rowV';
limit1 = 3;
limit2 = 1;
% bounds(1:3,1) = 1;
% bounds(1,2) = size(roiStack,1);
% bounds(2,2) = size(roiStack,2);
% bounds(3,2) = size(roiStack,3);
rlist = find(r(:,1)); %subtract from this list to remove used points
rowN = 1;
for j = 1:size(r,1) % run through every object
    if ismember(j,rlist) %check if current object has already been assign
        clear currentRow
        r(j,8) = rowN; %store row number
        currentRow = j; % set as first member in row
        rlist((j==rlist)) = 0; % remove j from list
        rLU = rlist;
        rLU(rLU==0) = max(rLU);
        rLU = unique(rLU);
        clear differences dv22
        
        %narrow potential row members with limit 1
        for k=1:size(rLU,1)
            if rLU(k,1) >0
                differences(k,1) = norm(cross(rowV,r(rLU(k,1),1:3)-r(j,1:3)))/norm(rowV);
            else
                differences(k,1)= 1000;
            end
        end
        rLUCand = rLU(differences(:,1)<limit1);
        
        %select first matches
        clear differences
        count = 0;
        for l = 1:size(rLUCand,1)
            if (rLUCand(l,1) >0) && rLUCand(l,1) ~= j
                differences(l,1) = norm(cross(rowV,r(rLUCand(l,1),1:3)-r(j,1:3)))/norm(rowV);
                count = count +1;
            else
                differences(l,1)= 1000;
                count = count +1;
            end
        end
        if count > 0
            rLUMatch = rLUCand(differences(:,1)<limit2);
            currentRow = cat(1,currentRow,rLUMatch);
        
        
        %loop until no additional matches are found
        
        matches = 1;
        while matches == 1
            matchFound = 0;
            for k = 1:size(rLUCand,1) %loop to select point of interest
                if rLUCand(k,1) > 0 && ismember(rLUCand(k,1),currentRow)==1
                    clear differences rLUMatch
                    for l = 1:size(rLUCand,1) %loop to create distances
                        if rLUCand(l,1) >0  && (rLUCand(l,1) ~= rLUCand(k,1))
                            differences(l,1) = norm(cross(rowV,r(rLUCand(l,1),1:3)-r(rLUCand(k,1),1:3)))/norm(rowV);
                        else
                            differences(l,1)= 1000;
                        end
                    end
                    rLUMatch = rLUCand(differences(:,1)<limit2);                   
                    rLUCand(k,1) = 0;
                    currentRow = cat(1,currentRow,rLUMatch);
                    if size(rLUMatch,1)>0
                        matchFound = 1;
                    end
                end
            end
            if matchFound == 0
                matches = 0;
            end
        end
        end
        
        
        
        currentRow = unique(currentRow);
        r(currentRow(currentRow>0),8) = rowN; %store row number
        rlist(currentRow(currentRow>0)) = 0; % remove j from list
        rows(rowN,1:size(currentRow,1)) = currentRow(:,1);
        rowN = rowN +1;
    end
end
end