function [r,rows] = buildRows2(r,rowV,planesFinal)

%scale rowV to length 2.12 in XY dimensions
rowV = rowV';
distX = 2.12; %microns
rowScale = sqrt(((distX)^2)/(rowV(1,1)^2+rowV(2,1)^2));
rowV = rowV*rowScale;

%set limits for row filters (microns)
limit1 = 4; %distance from line (rowV passing through j)
limit2 = 2; %distance from point (rowV+j)


rlist = find(r.X(:)); %subtract from this list to remove used points
rowN = 0;
for j = 1:r.l % run through every object
    if ismember(j,rlist) %check if current object has already been assign
        rowN = rowN+1;
        [jrow,jplane] = find(planesFinal==j);
        rlistP = planesFinal(1:nnz(planesFinal(:,jplane)),jplane);
        clear currentRow
        r.row(j) = rowN; %store row number
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
                differences(k,1) = norm(cross(rowV,r.r(rLU(k,1),1:3)-r.r(j,1:3)))/norm(rowV);
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
            differences(:,1) = sqrt(((r.X(j)*ones(size(rLUCand,1),1))+rowV(1,1)-r.X(rLUCand(:,1))).^2+((r.Y(j)*ones(size(rLUCand,1),1))+rowV(2,1)-r.Y(rLUCand(:,1))).^2);
            d1check = find(differences<limit2);
            if size(d1check,1)>0
                d1proceed = 1;
                seed = j;
                rowVS = 1;
                while d1proceed == 1
                    clear differences
                    differences(:,1) = sqrt(((r.X(seed)*ones(size(rLUCand,1),1))+(rowV(1,1)*rowVS)-r.X(rLUCand(:,1))).^2+(((r.Y(seed)*ones(size(rLUCand,1),1))+(rowV(2,1)*rowVS))-r.Y(rLUCand(:,1))).^2);
                    rLUMatch = find(differences<limit2);
                    if size(rLUMatch,1)>0
                        seed = rLUCand(find(differences==min(differences)));
                        r.row(seed) = rowN;
                        rlist(rlist(:,1)==seed,:) = [];
                        rowVS = 1;
                    else
                        rowVS =rowVS+1;
                        if r.X(seed)+(rowV(1,1)*rowVS)>max(r.X(:))+distX || r.Y(seed)+(rowV(2,1)*rowVS)>max(r.Y(:))+distX || r.X(seed)+(rowV(1,1)*rowVS)>min(r.X(:))-distX || r.Y(seed)+(rowV(2,1)*rowVS)>min(r.Y(:))-distX
                            d1proceed=0;
                        end
                    end
                    
                end
            end
            %direction2
            clear differences
            differences(:,1) = sqrt(((r.X(j)*ones(size(rLUCand,1),1))-rowV(1,1)-r.X(rLUCand(:,1))).^2+((r.Y(j)*ones(size(rLUCand,1),1))-rowV(2,1)-r.Y(rLUCand(:,1))).^2);
            d2check = find(differences<limit2);
            if size(d2check,1)>0
                d2proceed = 1;
                seed = j;
                rowVS = 1;
                while d2proceed == 1
                    clear differences
                    differences(:,1) = sqrt(((r.X(seed,1)*ones(size(rLUCand,1),1))-(rowV(1,1)*rowVS)-r.X(rLUCand(:,1),1)).^2+((r.Y(seed,1)*ones(size(rLUCand,1),1))-(rowV(2,1)*rowVS)-r.Y(rLUCand(:,1),1)).^2);
                    rLUMatch = find(differences<limit2);
                    if size(rLUMatch,1)>0
                        seed = rLUCand(find(differences==min(differences)));
                        r.row(seed) = rowN;
                        rlist(rlist(:,1)==seed,:) = [];
                        rowVS = 1;
                    else
                        rowVS =rowVS+1;
                        if r.X(seed)-(rowV(1,1)*rowVS)>max(r.X(:))+distX || r.Y(seed)-(rowV(2,1)*rowVS)>max(r.Y(:))+distX || r.X(seed)-(rowV(1,1)*rowVS)>min(r.X(:))-distX || r.Y(seed)-(rowV(2,1)*rowVS)>min(r.Y(:))-distX
                            d2proceed=0;
                        end
                    end                                        
                end
            end
        end
    end
end

for i = 1:rowN
    rows(i,1:size(find(r.row(:)==i))) = find(r.row(:)==i);
end