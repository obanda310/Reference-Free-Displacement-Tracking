function [r,rows] = buildRows(r,rowV,planesFinal)
rowV = rowV';


limit1 = 2.5;
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
        rLU(rLU==0,:) = [];
        if size(rLU,1) == 0
            clear rLU
            rLU = 0;
        end
        %rLU = unique(rLU);
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

%%
% maxX = max(r(:,1));
% maxY = max(r(:,2));
% if abs(rowV(1,1)) >abs(rowV(1,2))
%     wallStart = [1 2]; %X
% else
%     wallStart = [2 1]; %Y
% end
% [r,rOrig] =sortrows(r,wallStart);
% t = sqrt((2.18^2)/(rowV(1,1)^2+rowV(1,2)^2));
% if rowV(1,wallStart(1,1))<0
%     rowV = rowV*-t;
% else
%     rowV = rowV*t;
% end
% 
% r(:,8) = 0;
% rowN = 1;
% radXY = 1;
% 
% for i = 1:size(planesFinal,2)
%     for j = 1:size(r,1)
%         if r(j,8) == 0
%             r(j,8) = rowN;
%             [xCurr,yCurr] = find(planesFinal==j);
%             working = 1;
%             currLoc =  r(j,1:2) + rowV(1,1:2);
%             while working == 1
%                 clear distances rNbor3Val rNbor rNbor2 rNbor3
%                 botX = currLoc(1,1) -radXY;
%                 topX = currLoc(1,1) +radXY;
%                 botY = currLoc(1,2) -radXY;
%                 topY = currLoc(1,2) +radXY;
%                 if size(find(r(:,1)<topX & r(:,1)>botX & r(:,2)<topY & r(:,2)>botY & r(:,8)==0)) > 0
%                     rNbor(1,1:size(find(r(:,1)<topX & r(:,1)>botX & r(:,2)<topY & r(:,2)>botY & r(:,8)==0))) = find(r(:,1)<topX & r(:,1)>botX & r(:,2)<topY & r(:,2)>botY& r(:,8)==0);
%                     if size(intersect(rOrig(rNbor(:,1),1),planesFinal(:,yCurr)))>0
%                         rNbor2 = intersect(rOrig(rNbor(:,1),1),planesFinal(:,yCurr));
%                         for k = 1:size(rNbor2,1)
%                             rNbor3(k,1) = find(rOrig(:,1)==rNbor2(k,1));
%                         end
%                         rNbor3Val(:,1:2) = r(rNbor3,1:2);
%                         distances(:,1:2) = (rNbor3Val(:,1:2) - currLoc).^2;
%                         distances(:,3) = sqrt(distances(:,1)+distances(:,2));
%                         best = find(distances(:,3) == min(distances(:,3)));
%                         r(rNbor3(best,1),8) = rowN;
%                         currLoc =  r(rNbor3(best,1),1:2) + rowV(1,1:2);
%                     else
%                         currLoc = currLoc + rowV(1,1:2);
%                     end               
%                 else
%                         currLoc = currLoc + rowV(1,1:2);
%                 end
%                 if currLoc(1,1)>maxX+2 || currLoc(1,1)<-2 || currLoc(1,2)<-2 || currLoc(1,2)>maxY+2
%                     working = 0;
%                 end
%             end
%             rowN = rowN+1;
%         end
%     end
% end
% 
% [r4 orig] = sortrows(rOrig);
% r = r(orig,:);
% 
% for i=1:max(r(:,8))
%     num = size(find(r(:,8)==i),1);
% rows(i,1:num) = find(r(:,8)==i);
% end

end