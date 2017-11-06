function [r,rows] = buildRows2(r,rowV,planesFinal)

%scale rowV to length 2.12 in XY dimensions
rowV = rowV';
distX = 2.12; %microns
rowScale = sqrt(((distX)^2)/(rowV(1,1)^2+rowV(2,1)^2));
rowV = rowV*rowScale;

%set limits for row filters (microns)
limit1 = 4; %distance from line (rowV passing through j)
limit2 = 2; %distance from point (rowV+j)


rlist = find(r(:,1)); %subtract from this list to remove used points
rowN = 0;
for j = 1:size(r,1) % run through every object
    if ismember(j,rlist) %check if current object has already been assign
        rowN = rowN+1;
        [jrow,jplane] = find(planesFinal==j);
        rlistP = planesFinal(1:nnz(planesFinal(:,jplane)),jplane);
        clear currentRow
        r(j,8) = rowN; %store row number
        currentRow = j; % set as first member in row
        rlist((j==rlist)) = 0; % remove j from list
        clear rLU
        rLU = intersect(rlist,rlistP);
        rLU(rLU==0,:) = [];
        if size(rLU,1) == 0
            clear rLU
            rLU = 0;
        end
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
        
        
        %build rows feature by feature, assumption is that the closest
        %feature to j+/-rowV is the next row member
                
        if size(rLUCand,1)>0
            %direction1
            clear differences
            differences(:,1) = sqrt(((r(j,1)*ones(size(rLUCand,1),1))+rowV(1,1)-r(rLUCand(:,1),1)).^2+((r(j,2)*ones(size(rLUCand,1),1))+rowV(2,1)-r(rLUCand(:,1),2)).^2);
            d1check = find(differences<limit2);
            if size(d1check,1)>0
                d1proceed = 1;
                seed = j;
                rowVS = 1;
                while d1proceed == 1
                    clear differences
                    differences(:,1) = sqrt(((r(seed,1)*ones(size(rLUCand,1),1))+(rowV(1,1)*rowVS)-r(rLUCand(:,1),1)).^2+(((r(seed,2)*ones(size(rLUCand,1),1))+(rowV(2,1)*rowVS))-r(rLUCand(:,1),2)).^2);
                    rLUMatch = find(differences<limit2);
                    if size(rLUMatch,1)>0
                        seed = rLUCand(find(differences==min(differences)));
                        r(seed,8) = rowN;
                        rlist(rlist(:,1)==seed,:) = [];
                        rowVS = 1;
                    else
                        rowVS =rowVS+1;
                        if r(seed,1)+(rowV(1,1)*rowVS)>max(r(:,1))+distX || r(seed,2)+(rowV(2,1)*rowVS)>max(r(:,2))+distX || r(seed,1)+(rowV(1,1)*rowVS)>min(r(:,1))-distX || r(seed,2)+(rowV(2,1)*rowVS)>min(r(:,2))-distX
                            d1proceed=0;
                        end
                    end
                    
                end
            end
            %direction2
            clear differences
            differences(:,1) = sqrt(((r(j,1)*ones(size(rLUCand,1),1))-rowV(1,1)-r(rLUCand(:,1),1)).^2+((r(j,2)*ones(size(rLUCand,1),1))-rowV(2,1)-r(rLUCand(:,1),2)).^2);
            d2check = find(differences<limit2);
            if size(d2check,1)>0
                d2proceed = 1;
                seed = j;
                rowVS = 1;
                while d2proceed == 1
                    clear differences
                    differences(:,1) = sqrt(((r(seed,1)*ones(size(rLUCand,1),1))-(rowV(1,1)*rowVS)-r(rLUCand(:,1),1)).^2+((r(seed,2)*ones(size(rLUCand,1),1))-(rowV(2,1)*rowVS)-r(rLUCand(:,1),2)).^2);
                    rLUMatch = find(differences<limit2);
                    if size(rLUMatch,1)>0
                        seed = rLUCand(find(differences==min(differences)));
                        r(seed,8) = rowN;
                        rlist(rlist(:,1)==seed,:) = [];
                        rowVS = 1;
                    else
                        rowVS =rowVS+1;
                        if r(seed,1)-(rowV(1,1)*rowVS)>max(r(:,1))+distX || r(seed,2)-(rowV(2,1)*rowVS)>max(r(:,2))+distX || r(seed,1)-(rowV(1,1)*rowVS)>min(r(:,1))-distX || r(seed,2)-(rowV(2,1)*rowVS)>min(r(:,2))-distX
                            d2proceed=0;
                        end
                    end                                        
                end
            end
        end
    end
end

for i = 1:rowN
    rows(i,1:size(find(r(:,8)==i))) = find(r(:,8)==i);
end