classdef RowData
    properties
        m %members
        p %plane
        NDm %non-deformed members
        V %approximate vector
        VL %Vector Length (set to spacing on pattern,*need to update this to be adaptable)
        s %# of members
        Idx %a format which should become obsolete soon
        
    end
    methods
        %%
        function obj = RowData(r,raw,image) %initialize
            
            % This section of code uses user input to establish the row vector, which
            % should describe the most precise laser patterning direction. This vector
            % is used to build rows later.
            
            % First check to see if a previous row vector exists
            clear files check
            files = dir('Matlab Data Files\*.mat'); %Check Directory for default filenames
            if size(files,1)>=1
                for k = 1:size(files,1)
                    current=files(k).name;
                    check(k)=strcmp(current(end-7:end),'rowV.mat');
                end
                loc=find(check);
            else
                loc= zeros(1,1);
                loc(:,1) = [];
            end
            disp('Generating Row Slope')
            image = FindNDSquare(image);
            if size(loc,2)>0
                load('Matlab Data Files\rowV.mat')
                disp('Found a previous row slope!')
            else
                try
                [rowV] = rowVoptimizer(r,raw.pSpaceXY);
                save('Matlab Data Files\rowV.mat','rowV')   
                catch                
                [rowV] = rowVselector(image);
                save('Matlab Data Files\rowV.mat','rowV')
                end
                %rowV = rowV2;
            end
            rowV(1,3) = 0;
            obj.V = rowV;
            obj.V = obj.V';
            obj.VL = raw.pSpaceXY;%microns
            %scale rowV to length pSpaceXY in XY dimensions
            rowScale = sqrt(((obj.VL)^2)/(obj.V(1,1)^2+obj.V(2,1)^2));
            obj.V = obj.V*rowScale;
            disp(['done Generating Row Slope at ' num2str(toc) ' seconds'])
        end
        %%
        function [obj,r] = buildRows(obj,r,plane,raw,preview)
            disp('Building Rows')
            r.row = zeros(size(r.X));
            obj.m = [];
            rowN = 0;
            rlist = r.ND; %subtract from this list to remove used points
            for j = 1:r.l % run through every object
                if ismember(j,rlist) && r.row(j) == 0 %check if current object has been assigned
                    %set limits for row filters (microns)
                    limit1 = raw.pSpaceXY/3; %distance from line (rowV passing through j)
                    rowN = rowN+1;
                    [~,jplane] = find(plane.final==j);
                    jplane = unique(jplane);
                    rlistP = plane.final(1:nnz(plane.final(:,jplane)),jplane);
                    r.row(j) = rowN; %store row number
                    rlist(rlist==j,:) = []; % remove j from list
                    rLU = intersect(rlist,rlistP); %grab all non-assigned objects in the same plane that are not deformed
                    
                    %narrow potential row members with limit 1
                    RowLineDist = zeros(size(rLU));
                    for k=1:size(rLU,1)
                        tempd = zeros(1,3);
                        tempd(1,1:2) = r.r(rLU(k,1),1:2)-r.r(j,1:2);
                        RowLineDist(k,1) = norm(cross(obj.V,tempd))/norm(obj.V);
                    end
                    rLUCand = rLU(RowLineDist(:,1)<limit1); %final list of potential matches
                    r.row(rLUCand) = rowN;
                    for i = 1:size(rLUCand,1)
                        rlist(rlist==rLUCand(i,1),:) = [];
                    end
                end
            end
            
            %Second pass: Use ND objects to help grow rows from unassigned
            %objects. Use strict rules to prevent false matches
            rlist = r.D;
            limit1 = raw.pSpaceXY; %distance from main row line (rowV passing through Non-Deformed objects)
            limit2 = raw.pSpaceXY/5; %distance from point (rowV+j)
            for i = 1:max(r.row) %for each row,
                rowN = i;                
                iPlane = mode(r.plane(r.row==i));%determine plane
                
                if iPlane > 0 %ismember(iPlane,plane.groups(1,:)) == 0 %if plane is not in top plane group
                    rLU = intersect(rlist,plane.final(1:nnz(plane.final(:,iPlane)),iPlane)); %grab unmatched objects in that plane
                    
                    %determine which of those objects fall within the limit distance from rowN line
                    RowLineDist = zeros(size(rLU));
                    for k=1:size(rLU,1)
                        tempd = zeros(1,3);
                        tempd(1,1:2) = r.r(rLU(k,1),1:2)-mean(r.r((r.row==i),1:2),1); 
                        RowLineDist(k,1) = norm(cross(obj.V,tempd))/norm(obj.V); % distance of current unmatched object to line with slope 'rowV' passing through the average ND object position in current row
                    end
                    rLUCand = rLU(RowLineDist(:,1)<limit1); %final list of potential matches. i.e. objects within distance 'limit1' from line.
                    
                    %while there are still objects in this list, perform a
                    %point search and grab the next point which is closest to any
                    %point on the current row, and determine if its trajectory
                    %from the row line matches the closest assigned object in
                    %that row
                    if size(rLUCand,1) >0
                        proceed = 1;
                    else
                        proceed = 0;
                    end
                    while proceed == 1
                        rowList = find(r.row==i);
                        [rLIdx,rLDist] = dsearchn(r.r(rowList,1:2),r.r(rLUCand,1:2)); % find closest points in row to all points in potential candidates                        
          
                        thisCand = (rLUCand(rLDist==min(rLDist)));
                        thisIdx = rowList(rLIdx(rLDist==min(rLDist)));
                        thisDist = rLDist(rLDist==min(rLDist));
                        
                        if size(thisCand,1)>1 %in case of a tie of closest objects, pick the first one.
                            thisCand(2:end,:) = [];
                            thisIdx(2:end,:) = [];
                            thisDist(2:end,:) = [];
                        end
                        
                        tempd = zeros(1,3);
             
                        tempd(1,1:2) = r.r(thisCand,1:2)-mean(r.r((r.row==i)&r.NDl,1:2));
                        RLDist = norm(cross(obj.V,tempd))/norm(obj.V); %distance from unmatched to main row line
                        
                        tempd = zeros(1,3);
                        tempd(1,1:2) = r.r(thisCand,1:2)-(r.r(thisIdx,1:2));
                        RPDist = norm(cross(obj.V,tempd))/norm(obj.V); %distance from each other
                        
                        tempd = zeros(1,3);
                        tempd(1,1:2) = r.r(thisIdx,1:2)-mean(r.r((r.row==i)&r.NDl,1:2));
                        PPDist = norm(cross(obj.V,tempd))/norm(obj.V); %distance from matched to main row line
                        
                        
                        if ((RLDist < limit2 && PPDist < limit2) || (RPDist < limit2 && thisDist<(raw.pSpaceXY*1.5))) %if both are close to the line or to each other, assign this object to row
                            r.row(thisCand) = i;
                            rLUCand(rLUCand==thisCand,:) = [];%if trajectory is a match, assign it to row. Otherwise,
                        else
                            rLUCand(rLUCand==thisCand,:) = [];%remove it from this list and loop through again.
                        end
                        
                        
                        if size(rLUCand,1) ==0
                            proceed = 0;%if list is empty, break loop and transition to the next
                            %object
                        end
                        
                    end
                end
            end
            
            for i = 1:rowN
                obj.m(i,1:size(find(r.row(:)==i))) = find(r.row(:)==i);
            end
            
            if preview > 0
                figure
                hold on
                for i = 1:size(obj.m,1)
                    scatter3(r.X(obj.m(i,1:nnz(obj.m(i,:)))),r.Y(obj.m(i,1:nnz(obj.m(i,:)))),r.Z(obj.m(i,1:nnz(obj.m(i,:)))))
                end
            end
        end
        
        %%
        function [obj,r]=formatRows(obj,plane,r)
            obj.p = [];
            % Separate Rows by Plane
            clear rowPlanes
            for i = 1:size(plane.groups,1)
                for j = 1:size(obj.m,1)
                    tIdx = [];
                    pIdx = plane.groups(i,1:nnz(plane.groups(i,:)));
                    for k = 1:nnz(pIdx)
                    tIdx = cat(1,tIdx,plane.final(:,pIdx(k)));
                    end
                    tIdx(tIdx(:,1)==0,:) = [];
                    rowPlanes(j,1:nnz(intersect(obj.m(j,:),tIdx)),i) = intersect(obj.m(j,:),tIdx);
                end
            end
            
            % Remake 'rows' variable with rows ordered by plane
            clear newRows
            nRS = find(rowPlanes(:,1,1)>0);
            newRows(:,:) = rowPlanes(nRS,:,1);
            rowPlanesIdx(1,1) = 1;
            rowPlanesIdx(1,2) = size(nRS,1);
            for i = 2:size(rowPlanes,3)
                clear currentPlanes
                nRS = find(rowPlanes(:,1,i)>0,1,'first');
                currentPlanes(:,:) =  rowPlanes((rowPlanes(:,1,i)>0),:,i);
                rowPlanesIdx(i,1) = size(newRows,1)+1;
                newRows = cat(1,newRows,currentPlanes);
                rowPlanesIdx(i,2) = size(newRows,1);
            end
            obj.m = newRows;
            obj.Idx = rowPlanesIdx;
            
            for i = 1:size(rowPlanesIdx,1)
                obj.p(rowPlanesIdx(i,1):rowPlanesIdx(i,2),1) = i;
            end
            
            % Label Objects in 'r' with Their Row Number
            for i = 1:size(obj.m,1)
                for j = 1:size(obj.m,2)
                    if obj.m(i,j)>0
                        r.row(obj.m(i,j),1) = i;
                    end
                end
            end
            
            % Determine non-deformed members of each row
            clear rowsNDCU rowsNDC
            rowsNDC = obj.m;
            for i = 1:size(obj.m,1)
                for j = 1:size(obj.m,2)
                    if ismember(obj.m(i,j),r.ND) == 0
                        rowsNDC(i,j) = 0;
                    end
                end
                rowsNDC(i,find(rowsNDC(i,:)==0)) = max(rowsNDC(i,:));
                rowsNDCU(i,1:size(unique(rowsNDC(i,:)),2)) = unique(rowsNDC(i,:));
            end
            obj.NDm = rowsNDCU;
            
            %Sort rows by X and Y values using row vector
            for i =1:size(obj.m,1)
                temprow = zeros(size(obj.m(i,1:nnz(obj.m(i,:)))));
                temprow(1,:) = obj.m(i,1:nnz(obj.m(i,:)));
                tlocs = zeros(size(temprow,2),2);
                tlocs(:,1) = r.X(temprow);
                tlocs(:,2) = r.Y(temprow);
                if abs(obj.V(1,1))<abs(obj.V(2,1))
                    [~,order] = sortrows(tlocs,[2,1]);
                else
                    [~,order] = sortrows(tlocs,[1,2]);
                end
                obj.m(i,1:nnz(obj.m(i,:))) = temprow(1,order);
            end
            
            
        end
        function obj = rowSizes(obj)
            obj.s = [];
            for i = 1:size(obj.m,1)
                obj.s(i,1) = nnz(obj.m(i,:));
            end
        end
        function ViewRows(obj,r)
            figure
            hold on
            for i = 1:size(obj.m,1)
                plot3(r.X(obj.m(i,1:nnz(obj.m(i,:)))),r.Y(obj.m(i,1:nnz(obj.m(i,:)))),r.Z(obj.m(i,1:nnz(obj.m(i,:)))))
            end
        end
    end
end