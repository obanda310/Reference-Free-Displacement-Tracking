classdef VertPlaneData
    properties
        raw
        rows %Rows of objects in this V-Plane
        cols %Columns of objects in this V-Plane
        colsO %Columns of objects in this V-Plane sorted by XY position
        planes %Planes in this V-Plane
        pgs %plane group idx (see planedata properties)
        grid %indices of detections in gridded format (row,col,vertplane)
        gridRow
        gridCol
        griddS %shear displacement (same format as 'grid')
        griddX
        griddY
        vX % d1 = vplane, d2 = plane, d3 = xyz components
        vY
        vZ
        
    end
    methods
        function [obj,r] = VertPlaneData(r,rows)
            clear RC
            RC(:,1) = r.row(:,1);
            RC(:,2) = r.col(:,1);
            RC2 = sortrows(RC,[2 1]);
            RC2(RC2(:,2)==0,:) = [];
            RC2(RC2(:,1)==0,:) = [];
            
            clear RCPlanes
            for i = 1:max(RC2(:,2))
                temp = RC2(RC2(:,2)==i,1);
                temp(temp==0,:) = [];
                RCPlanes(i,1:nnz(temp)) = temp;
            end
            obj.raw = unique(RCPlanes,'rows'); %get rid of duplicate planes            
            obj = cleanRows(obj,r);            
            
            for i = 1:size(obj.raw,1)
                for j = 1:size(obj.raw,2)
                    if j == 1
                        colList = r.col(r.row == obj.raw(i,j));
                    elseif obj.raw(i,j)>0
                        colList = cat(1,colList,r.col(r.row == obj.raw(i,j)));
                    end
                end
                colListU = unique(colList); %these are all columns represented in current vplane
                colListU(colListU==0,:) = [];
                
                %now use the columns to count how many instances of each
                %row exist in the vplane. If there are rows that occur very
                %infrequently, identify which columns are introducing these
                %rows and remove the column information from the offending
                %objects
                for j = 1:size(colListU)
                    if j == 1
                        rowList = r.row(r.col == colListU(j));
                        colList2 = r.col(r.col == colListU(j));
                    else
                        rowList = cat(1,rowList,r.row(r.col == colListU(j)));
                        colList2 = cat(1,colList2,r.col(r.col == colListU(j)));
                    end
                end
                
                rowListU = unique(rowList);
                rowListU(rowListU==0,:) = [];
                
                
                for j = 1:nnz(rowListU)
                    rowListP(j) = mode(r.plane(r.row==rowListU(j)));
                end
                
                
                clear rowCounts
                for j = 1:nnz(rowListU)                    
                    if nnz(find(rowListP==rowListP(j)))>1 % before removing any column data, make sure there will still be at least one row in each plane!
                        rowCounts(j,1) = nnz(rowList == rowListU(j));
                        if rowCounts(j,1) < 5
                            rowListP(j) = 0;
                            colIdx = colList2(rowList == rowListU(j));
                            for k = 1:size(colIdx,1)
                                r.col(r.row==rowListU(j) & r.col == colIdx(k)) = 0;
                            end
                        end
                    end
                end
                
            end
            
            clear RC
            RC(:,1) = r.row(:,1);
            RC(:,2) = r.col(:,1);
            RC2 = sortrows(RC,[2 1]);
            RC2(RC2(:,2)==0,:) = [];
            
            clear RCPlanes
            for i = 1:max(RC2(:,2))
                temp = RC2(RC2(:,2)==i,1);
                temp(temp==0,:) = [];
                RCPlanes(i,1:nnz(temp)) = temp;
            end
            obj.raw = unique(RCPlanes,'rows'); %get rid of duplicate planes
            obj = cleanRows(obj,r);
            
            
            % Sort planes based on physical locations
            for i = 1:size(obj.raw,1)
                for j = 1:size(obj.raw(i,:),2)
                    if j == 1 && obj.raw(i,j) > 0
                        Xs = r.X(r.row==obj.raw(i,j));
                        Ys = r.Y(r.row==obj.raw(i,j));
                    else
                        Xs = cat(1,Xs,r.X(r.row==obj.raw(i,j)));
                        Ys = cat(1,Ys,r.Y(r.row==obj.raw(i,j)));
                    end
                end
                meanXY(i,1) = mean(Xs);
                meanXY(i,2) = mean(Ys);
            end
            if abs(rows.V(1,1))<abs(rows.V(2,1))
                [~,order] = sortrows(meanXY,[1,2]);
            else
                [~,order] = sortrows(meanXY,[2,1]);
            end
            obj.raw = obj.raw(order,:);
            obj.raw(isnan(obj.raw))= 0;
            
            
            function obj = cleanRows(obj,r)
                for i2 = 1:size(obj.raw,1)
                    clear temp2
                    temp2(:,1) = obj.raw(i2,:);
                    temp2 = unique(temp2);
                    temp2(isnan(temp2),:) = [];
                    temp2(temp2==0,:)=[];                    
                    if size(temp2,1)>0
                        obj.raw(i2,:) = 0;
                        obj.raw(i2,1:nnz(temp2)) = temp2;
                        tidx = [];
                        for j2 = 1:nnz(obj.raw(i2,:))
                            if obj.raw(i2,j2)>0
                                tidx = [tidx; find(r.row==obj.raw(i2,j2))];
                            end
                        end                       
                        pNum(i2) = max(r.colSize(tidx));
                        
                        if nnz(obj.raw(i2,:))<pNum(i2) % need to update this!!
                            obj.raw(i2,:) = 0;
                        end
                    end
                end
                obj.raw(isnan(obj.raw)) = 0;
                obj.raw(obj.raw(:,1)==0,:) = [];
                obj.raw = unique(obj.raw,'rows'); %get rid of duplicate planes
            end
            
        end
        
        
        function obj = removeDuplicates(obj)
            % Intersect each 'plane' with every other to determine if repeats exist
            vpOverlap = zeros(size(obj.raw,1));
            for i = 1:size(obj.raw,1)
                for j = 1:size(obj.raw,1)
                    if j ~= i
                        check = intersect(obj.raw(i,1:nnz(obj.raw(i,:))),obj.raw(j,1:nnz(obj.raw(j,:))));
                        if size(check,2)>0
                            vpOverlap(i,j) = 1;
                        end
                    end
                end
            end
            
            
            if sum(vpOverlap(:))>0
                vpOverlap = vpOverlap>0;
                vpIdx = 1:size(obj.raw,1);
                for i = 1:size(obj.raw,1)
                    clear overlapIdx vpCurrent
                    overlapIdx = vpIdx(vpOverlap(i,:));
                    vpCurrent = obj.raw(i,:)';
                    
                    if size(overlapIdx,2)>0
                        for j = 1:size(overlapIdx,2)
                            vpCurrent = cat(1,vpCurrent,obj.raw(overlapIdx(j),:)');
                        end
                        vpCurrent = unique(vpCurrent);
                    end
                    vpCurrent(vpCurrent==0,:) = []; %get rid of zeros
                    obj.rows(i,1:nnz(vpCurrent)) = vpCurrent;
                end
                obj.rows = unique(obj.rows,'rows');
            else
                obj.rows=obj.raw;
            end
        end
        
        function [obj,r,rows,eraseCheck] = cleanVPlanes(obj,r,rows)
            % Identify groups with multiple rows per plane and merge co-planar rows
            % Then update r and rows to reflect merged rows
            eraseCheck = 1;
            for i = 1:size(obj.rows,1)
                for j = 1:size(obj.rows,2)
                    if obj.rows(i,j) > 0
                        RCPF2(i,j) = rows.p(obj.rows(i,j));
                    end
                end
            end
            %
            for i = 1:size(RCPF2)
                if length(RCPF2(i,1:nnz(RCPF2(i,:)))) ~= length(unique(RCPF2(i,1:nnz(RCPF2(i,:)))))
                    RCPFchk(i,1) = 1;
                else
                    RCPFchk(i,1) = 0;
                end
            end
            RCPFchk2 = find(RCPFchk);
            count = 1;
            clear erasedRows
            erasedRows = 0;
            for i2 = 1:size(RCPFchk2)
                i = RCPFchk2(i2);
                for j = 2:nnz(RCPF2(i,:))
                    if RCPF2(i,j) == RCPF2(i,j-1)
                        EndIdx = nnz(rows.m(obj.rows(i,j-1),:));
                        rows.m(obj.rows(i,j-1),EndIdx+1:EndIdx+nnz(rows.m(obj.rows(i,j),:))) = rows.m(obj.rows(i,j),1:nnz(rows.m(obj.rows(i,j),:)));
                        for k = 1:nnz(rows.m(obj.rows(i,j),:))
                            r.row(rows.m(obj.rows(i,j),k)) = obj.rows(i,j-1);
                        end
                        rows.m(obj.rows(i,j),1:nnz(rows.m(obj.rows(i,j),:))) = zeros(size(1:nnz(rows.m(obj.rows(i,j),:))));
                        %record which rows are being moved before erasing
                        %them
                        erasedRows(count) = obj.rows(i,j);
                        count = count+1;
                        obj.rows(i,j) = obj.rows(i,j-1);
                    end
                end
            end
            % Next, clean all locations that referenced old row structure
            % This includes: rows.m,r.rows,
            rows.m(rows.m(:,1)==0,:) = [];
            if size(erasedRows,1) == 1 && erasedRows(1,1) == 0
                eraseCheck = 0;
            else
                for i = 1:size(erasedRows)
                    r.row(r.row(:,1)>erasedRows(i),1) = r.row(r.row(:,1)>erasedRows(i),1)-1;
                end
            end
            
        end
        
        
        function [obj,r] = ProfileColumns(obj,r,plane)
            %% Determine Average Pillar in VPlanes (these will be used to calculate ref. positions.)
            % Should identify the average vectors to get from one non-deformed object
            % within a CR Plane to the one above it in the same CR Plane.
            
            %for each vplane
            obj.cols = [];
            obj.planes = [];
            
            r.colCheck = zeros(size(r.X));
            r.colUp = zeros(size(r.X));
            
            
            for i = 1:max(r.vplane)
                
                vpp = unique(r.plane(r.vplane==i));%determine which planes are in vplane
                vpp(vpp==0,:) = [];
                
                %determine the order of the planes
                vpplocs = zeros(size(vpp));
                for j = 1:size(vpp)
                    vpplocs(j,1) = plane.loc(vpp(j,1));
                end
                
                [vpplocs, idx]=sortrows(vpplocs);
                vpp = vpp(idx); %sort the planes in this vplane by Z location
                
                obj.planes(i,1:nnz(vpp)) = vpp; %store the sorted planes                
                %find all members of current vplane and determine if the next member in
                %their column is in the next plane
                temp = unique(r.col(r.vplane(:,1)==i));
                temp(temp==0,:) = [];
                
                obj.cols(i,1:nnz(temp)) = temp; % store this vplane's columns
                for j = 1:nnz(obj.cols(i,:)) %for each column,
                    thisColPlanes = r.plane(r.col(:,1)==obj.cols(i,j),1); % find all planes in the column
                    for k = 1:size(vpp,1)-1
                        if ismember(vpp(k,1),thisColPlanes) && ismember(vpp(k+1,1),thisColPlanes)
                            r.colUp(r.col(:,1)==obj.cols(i,j) & r.plane(:,1)==vpp(k,1),1) = find(r.col(:,1)==obj.cols(i,j) & r.plane(:,1)==vpp(k+1,1));
                            if min(r.NDl(r.col(:,1)==obj.cols(i,j) & r.plane(:,1)==vpp(k,1),1)) > 0 &&  min(r.NDl(r.col(:,1)==obj.cols(i,j) & r.plane(:,1)==vpp(k+1,1),1)) > 0 %make sure they are in non-deformed region
                                r.colCheck(r.col(:,1)==obj.cols(i,j) & r.plane(:,1)==vpp(k,1),1) = 1;
                                if k == size(vpp,1)-1
                                    r.colCheck(r.col(:,1)==obj.cols(i,j) & r.plane(:,1)==vpp(k+1,1),1) = 1;
                                end
                            end
                        end
                    end
                end
                %end
            end
            
            obj.vX = [];
            obj.vY = [];
            obj.vZ = [];
            for i = 1:max(r.vplane)
                %find average vectors for that vplane
                for j = 1:nnz(obj.planes(i,:))-1
                    obj.vX(i,j) = (sum(r.X(r.colUp(r.colCheck(:,1)==1 & r.vplane(:,1) == i & r.plane(:,1)==obj.planes(i,j),1),1) - r.X(r.colCheck(:,1)==1 & r.vplane(:,1) == i & r.plane(:,1)==obj.planes(i,j),1)))/nnz(r.X(r.colCheck(:,1)==1 & r.vplane(:,1) == i & r.plane(:,1)==obj.planes(i,j),1));
                    obj.vY(i,j) = (sum(r.Y(r.colUp(r.colCheck(:,1)==1 & r.vplane(:,1) == i & r.plane(:,1)==obj.planes(i,j),1),1) - r.Y(r.colCheck(:,1)==1 & r.vplane(:,1) == i & r.plane(:,1)==obj.planes(i,j),1)))/nnz(r.Y(r.colCheck(:,1)==1 & r.vplane(:,1) == i & r.plane(:,1)==obj.planes(i,j),1));
                    obj.vZ(i,j) = (sum(r.Z(r.colUp(r.colCheck(:,1)==1 & r.vplane(:,1) == i & r.plane(:,1)==obj.planes(i,j),1),1) - r.Z(r.colCheck(:,1)==1 & r.vplane(:,1) == i & r.plane(:,1)==obj.planes(i,j),1)))/nnz(r.Z(r.colCheck(:,1)==1 & r.vplane(:,1) == i & r.plane(:,1)==obj.planes(i,j),1));
                end
            end
            
            obj.pgs = zeros(size(obj.planes));
            for i = 1:size(obj.planes,1)
                for j = 1:size(obj.planes,2)
                    if obj.planes(i,j) >0
                        [rowIdx,colIdx] = find(plane.groups==obj.planes(i,j));
                        obj.pgs(i,j) = rowIdx;
                    end
                end
            end
        end
        
        
        function [obj,r] = AddReferenceLocs(obj,r)
            % these values are temporary and will be replaced once rows have been fit to a line
            r.rX = zeros(size(r.X));
            r.rY = zeros(size(r.X));
            r.rZ = zeros(size(r.X));
            for i = 1:size(r.X,1)
                if r.colUp(i) > 0
                    if r.vplane(i) > 0
                        r.rX(r.colUp(i),1) = obj.vX(r.vplane(i),obj.planes(r.vplane(i),:)==r.plane(i))+r.X(i);
                        r.rY(r.colUp(i),1) = obj.vY(r.vplane(i),obj.planes(r.vplane(i),:)==r.plane(i))+r.Y(i);
                        r.rZ(r.colUp(i),1) = obj.vZ(r.vplane(i),obj.planes(r.vplane(i),:)==r.plane(i))+r.Z(i);
                        
                    elseif r.vplane(i) == 0 && r.vplane(r.colUp(i))>0
                        r.vplane(i) = r.vplane(r.colUp(i));
                        r.rX(r.colUp(i),1) = obj.vX(r.vplane(i),obj.planes(r.vplane(i),:)==r.plane(i))+r.X(i);
                        r.rY(r.colUp(i),1) = obj.vY(r.vplane(i),obj.planes(r.vplane(i),:)==r.plane(i))+r.Y(i);
                        r.rZ(r.colUp(i),1) = obj.vZ(r.vplane(i),obj.planes(r.vplane(i),:)==r.plane(i))+r.Z(i);
                    end
                end
            end
            % The bottom of each pillar should now have r.colUp>0 and r.rX == 0.
            % Isolate those and set reference locations equal to actual locations
            
            % again, these values are temporary
            r.rX(r.colUp>0 & r.rX==0) = r.X(r.colUp>0 & r.rX==0);
            r.rY(r.colUp>0 & r.rY==0) = r.Y(r.colUp>0 & r.rY==0);
            r.rZ(r.colUp>0 & r.rZ==0) = r.Z(r.colUp>0 & r.rZ==0);
        end
        
        function [obj,r] = BuildGrids(obj,r,rows)
            %% Use catalogued data to grid available data and identify missing data
            % Build out each vertical plane in a gridded format (i.e. matrix rows =
            % actual rows, matrix columns = actual columns). Empty spaces (zeros)should
            % denote missing information.
            obj.colsO = zeros(size(obj.cols));
            obj.grid = zeros(max(obj.pgs(:)),size(rows.m,2),size(obj.cols,1));
            obj.gridRow = zeros(size(obj.grid));
            obj.gridCol = zeros(size(obj.grid));
            r.rowLeft =zeros(size(r.X));
            r.rowRight = zeros(size(r.X));
            r.colUp = zeros(size(r.X));
            r.colDown = zeros(size(r.X));
            
            for i = 1:size(obj.cols,1) %for each vplane
                % first sort each row by x then y
                clear tcols tlocs
                tcols = obj.cols(i,1:nnz(obj.cols(i,:)))';
                if size(tcols,1)>0
                    for j = 1:size(tcols,1)
                        tlocs(j,1) = r.X(r.col==tcols(j) & r.plane == max(r.plane(r.col==tcols(j)))); %mean x pos
                        tlocs(j,2) = r.Y(r.col==tcols(j) & r.plane == max(r.plane(r.col==tcols(j))));  %mean y pos
                    end
                    if abs(rows.V(1,1))<abs(rows.V(2,1))
                        [~,order] = sortrows(tlocs,[2,1]);
                    else
                        [~,order] = sortrows(tlocs,[1,2]);
                    end
                    tcols = tcols(order)';
                    obj.colsO(i,1:size(tcols,2),1) = tcols;
                    
                    % then use the sorted column to fill out the full grid.
                    for j = 1:size(tcols,2)
                        for k = 1:max(obj.pgs(:))
                            if size(find(r.planeGroup==k&r.col==tcols(1,j)),1)==1
                                obj.grid(k,j,i) = find(r.planeGroup==k&r.col==tcols(1,j));
                                obj.gridCol(k,j,i) = r.col(find(r.planeGroup==k&r.col==tcols(1,j)));
                                obj.gridRow(k,j,i) = r.row(find(r.planeGroup==k&r.col==tcols(1,j)));
                            end
                        end
                    end
                else
                end
            end
            % Replace all empty "columns" with NaN
            for i = 1:size(obj.grid,3)
                for j = 1:size(obj.grid,2)
                    if obj.grid(:,j,i) == zeros(size(obj.grid,1),1)
                        obj.grid(:,j,i) = NaN;
                    end
                end
            end
            % Update neighbor information for each object
            for i = 1:size(obj.grid,1)
                for j = 1:size(obj.grid,2)
                    for k = 1:size(obj.grid,3)
                        cidx = obj.grid(i,j,k);
                        if cidx>0
                            %left
                            if j>1
                                if isnan(obj.grid(i,j-1,k)) ==0
                                    r.rowLeft(cidx) = obj.grid(i,j-1,k);
                                end
                            end
                            %right
                            if j<size(obj.grid,2)
                                if isnan(obj.grid(i,j+1,k)) ==0
                                    r.rowRight(cidx) = obj.grid(i,j+1,k);
                                end
                            end
                            %bottom
                            if i<size(obj.grid,1)
                                r.colDown(cidx) = obj.grid(i+1,j,k);
                            end
                            %top
                            if i>1
                                r.colUp(cidx) = obj.grid(i-1,j,k);
                            end
                        end
                    end
                end
            end
        end
        %%
        function obj = BuildDispGrids(obj,r)
            % Build grids of displacement data
            obj.griddS = zeros(size(obj.grid));
            obj.griddX = zeros(size(obj.grid));
            obj.griddY = zeros(size(obj.grid));
            for i = 1:size(obj.grid,3)
                for j = 1:size(obj.grid,2)
                    for k =1:size(obj.grid,1)
                        if obj.grid(k,j,i) > 0 && isnan(obj.grid(k,j,i)) ==0 %&& ismember(obj.grid(k,j,i),r.D)==1 %&& m3.disp(obj.grid(k,j,i),4)>m3.SnoiseCO
                            obj.griddS(k,j,i) = r.dS(obj.grid(k,j,i));
                        end
                        if obj.grid(k,j,i) > 0 && isnan(obj.grid(k,j,i)) ==0 %&& ismember(obj.grid(k,j,i),r.D)==1 %&& m3.disp(obj.grid(k,j,i),1)>m3.XnoiseCO
                            obj.griddX(k,j,i) = r.dX(obj.grid(k,j,i));
                        end
                        if obj.grid(k,j,i) > 0 && isnan(obj.grid(k,j,i)) ==0 %&& ismember(obj.grid(k,j,i),r.D)==1 %&& m3.disp(obj.grid(k,j,i),2)>m3.YnoiseCO
                            obj.griddY(k,j,i) = r.dY(obj.grid(k,j,i));
                        end
                    end
                end
            end
            
            for i = 1:size(obj.griddX,3)
                for j = 1:size(obj.griddX,2)
                    if obj.griddX(:,j,i) == zeros(size(obj.grid,1),1)
                        obj.griddX(:,j,i) = NaN;
                        obj.griddY(:,j,i) = NaN;
                        obj.griddS(:,j,i) = NaN;
                    end
                end
            end
        end
        %%
        %Use row information to fill missing column slots
        function [obj,r] = FillHoles(obj,r,limit)
            
            %col holes
            tidx1 = find(r.row>0 & r.col == 0); %objects to fill holes
            tidx2 = find(r.colUp == 0 & r.planeGroup ~=1); %holes to be filled
            if nnz(tidx1) >0 && nnz(tidx2)>0
                [pIdx,pDist] = dsearchn(r.r(tidx1,1:2),r.r(tidx2,1:2));
                for i = 1:nnz(tidx2)
                    if pDist(i)< limit %&& (r.vplane(tidx1(pIdx)) == r.vplane(tidx2(i)))
                        targetRow = r.row(tidx1(pIdx(i)));
                        [targetVP,~] = find(obj.rows==targetRow);
                        thisRow = r.row(tidx2(i));
                        [thisVP,~] = find(obj.rows==thisRow);
                        if targetVP == thisVP
                            r.col(tidx1(pIdx(i))) = r.col(tidx2(i));
                            r.vplane(tidx1(pIdx(i))) = r.vplane(tidx2(i));
                        end
                    end
                end
                
                tidx1 = find(r.col>0 & r.row == 0 & r.vplane>0); %objects to fill holes
                for i = 1:nnz(tidx1)
                    targetPlane = r.planeGroup(tidx1(i));
                    targetVP = r.vplane(tidx1(i));
                    targetRow = mode(obj.gridRow(targetPlane,:,targetVP));
                    r.row(tidx1(i)) = targetRow;
                end
            end
        end
        %% Reset certain objects in V-Plane
        function obj = resetGridObject(obj,idx)
            for i = 1:size(idx,1)
                obj.gridRow(obj.grid==(idx(i))) = 0;
                obj.gridCol(obj.grid==(idx(i))) = 0;
                obj.griddS(obj.grid==(idx(i))) = 0;
                obj.griddX(obj.grid==(idx(i))) = 0;
                obj.griddY(obj.grid==(idx(i))) = 0;
                obj.grid(obj.grid==(idx(i))) = 0;
                
            end
        end
            %%
            %for diagnostics
        function ViewVPlanes(obj,r)
            tempmap = brewermap(size(obj.grid,3)*2,'GnBu');
            figure
            hold on
            for i = 1:size(obj.grid,3)
                temp = obj.grid(:,:,i);
                temp2 = temp(:);
                temp2(temp2==0,:) = [];
                temp2(isnan(temp2),:) = [];
                plot3(r.X(temp2),r.Y(temp2),r.Z(temp2),'Color',tempmap(i+round(.25*size(obj.grid,3)),:))
            end
            
            tIdx = r.row==0 | r.col==0;
            if sum(tIdx) >0
                scatter3(r.X(tIdx),r.Y(tIdx),r.Z(tIdx),20,'red','x')
            end
            
            tIdx = r.colUp==0 & r.planeGroup~=1;
            if sum(tIdx) >0
                scatter3(r.X(tIdx),r.Y(tIdx),r.Z(tIdx),20,'blue','*')
            end            
            
            if r.State>1
                quiver3(r.rX,r.rY,r.rZ,r.dX,r.dY,r.dZ,0,'Color','blue')
            end
        end
        
        %% for publication
        function ViewVPlanes2(obj,r)
            tempmap = brewermap(size(obj.grid,3)*2,'GnBu');
            vp=figure;
            hold on
            for i = 1:size(obj.grid,3)
                temp = obj.grid(:,:,i);
                temp2 = temp(:);
                temp2(temp2==0,:) = [];
                temp2(isnan(temp2),:) = [];
                scatter3(r.X(temp2),r.Y(temp2),r.Z(temp2),'.') %,'Color',tempmap(i+round(.25*size(obj.grid,3)),:)
            end
            
%             tIdx = r.row==0 | r.col==0;
%             if sum(tIdx) >0
%                 scatter3(r.X(tIdx),r.Y(tIdx),r.Z(tIdx),20,'red','x')
%             end
%             
%             tIdx = r.colUp==0 & r.planeGroup~=1;
%             if sum(tIdx) >0
%                 scatter3(r.X(tIdx),r.Y(tIdx),r.Z(tIdx),20,'blue','*')
%             end
%             
%             if r.State>1
%                 quiver3(r.rX,r.rY,r.rZ,r.dX,r.dY,r.dZ,0,'Color','blue')
%             end
            
%             bcolor = 'white';
%             fcolor = 'black';
            bcolor = 'black';
            fcolor = 'white';
            AxisFontSize = 12;
            LegendFontSize = 14;
            xt = 'Y \mum';% input('enter the xaxis label','s');
            yt = 'X \mum'; %input('enter the yaxis label','s');
            zt = 'Z \mum';
            label{1} = xlabel(xt);
            label{2} = ylabel(yt);
            label{3} = zlabel(zt);
            set(gca,'YMinorTick','on','color',bcolor)
            ytickformat('%.1f')
            le{1} = 'plane 1'; %input('enter the legend','s');
            le{2} = 'plane 2';
            ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,1,[0 0])
            %errorbar(meanDisplacements(1,1:3),meanDisplacements(2,1:3),'.','color',[0 0 0],'MarkerSize',1)
            axis([0 max(r.X) 0 max(r.Y) 0 ceil(max(r.Z))])
            legend off
            hold off
            view(-15,15)
            savefile = 'vplanes.tif';
            export_fig(vp,savefile,'-native');
        end
                %% for publication
        function ViewVPlanes3(obj,r)
            tempmap = brewermap(size(obj.grid,3)*2,'GnBu');
            vp=figure;
            hold on
            for i = 1:size(obj.grid,3)
                temp = obj.grid(:,:,i);
                temp2 = temp(:);
                temp2(temp2==0,:) = [];
                temp2(isnan(temp2),:) = [];
                scatter3(r.X(temp2),r.Y(temp2),r.Z(temp2),5,'.') %,'Color',tempmap(i+round(.25*size(obj.grid,3)),:)
                plot3(r.X(temp2),r.Y(temp2),r.Z(temp2),'Color',tempmap(i+round(.25*size(obj.grid,3)),:))
            end
            
            plot3([10 30 30 10 10],[25 25 60 60 25],[0 0 0 0 0],'g')
            
            tIdx = r.row==0 | r.col==0;
            if sum(tIdx) >0
                scatter3(r.X(tIdx),r.Y(tIdx),r.Z(tIdx),20,'green','x')
            end
            
            tIdx = r.colUp==0 & r.planeGroup~=1;
            if sum(tIdx) >0
                scatter3(r.X(tIdx),r.Y(tIdx),r.Z(tIdx),20,'cyan','*')
            end
            
            if r.State>1
                quiver3(r.rX,r.rY,r.rZ,r.dX,r.dY,r.dZ,0,'Color','blue')
            end
            
%             bcolor = 'white';
%             fcolor = 'black';
            bcolor = 'black';
            fcolor = 'white';
            AxisFontSize = 12;
            LegendFontSize = 14;
            xt = 'X \mum';% input('enter the xaxis label','s');
            yt = 'Y \mum'; %input('enter the yaxis label','s');
            zt = 'Z \mum';
            label{1} = xlabel(xt);
            label{2} = ylabel(yt);
            label{3} = zlabel(zt);
            set(gca,'YMinorTick','on','color',bcolor)
            ytickformat('%.1f')
            le{1} = 'plane 1'; %input('enter the legend','s');
            le{2} = 'plane 2';
            ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,1,[0 0])
            %errorbar(meanDisplacements(1,1:3),meanDisplacements(2,1:3),'.','color',[0 0 0],'MarkerSize',1)
            axis([0 max(r.X) 0 max(r.Y) 0 ceil(max(r.Z))])
            legend off
            hold off
            view(0,90)
            savefile = 'vplanes3.tif';
            export_fig(vp,savefile,'-native');
        end
                        %% for publication
        function ViewVPlanes4(obj,r)
            tempmap = brewermap(size(obj.grid,3)*2,'GnBu');
            vp=figure('Position', [10 10 900 900]);
            
            hold on
            for i = 1:size(obj.grid,3)
                temp = obj.grid(:,:,i);
                temp2 = temp(:);
                temp2(temp2==0,:) = [];
                temp2(isnan(temp2),:) = [];
                scatter3(r.X(temp2),r.Y(temp2),r.Z(temp2),'.') %,'Color',tempmap(i+round(.25*size(obj.grid,3)),:)
                %plot3(r.X(temp2),r.Y(temp2),r.Z(temp2),'Color',tempmap(i+round(.25*size(obj.grid,3)),:))
            end
            
            tIdx = r.row==0 | r.col==0;
            if sum(tIdx) >0
                scatter3(r.X(tIdx),r.Y(tIdx),r.Z(tIdx),20,'white','x')
            end
            
            tIdx = r.colUp==0 & r.planeGroup~=1;
            if sum(tIdx) >0
                scatter3(r.X(tIdx),r.Y(tIdx),r.Z(tIdx),20,'green','*')
            end
            
            if r.State>1
                quiver3(r.rX,r.rY,r.rZ,r.dX,r.dY,r.dZ,0,'Color','blue')
            end
            
%             bcolor = 'white';
%             fcolor = 'black';
            bcolor = 'black';
            fcolor = 'white';
            AxisFontSize = 12;
            LegendFontSize = 14;
            xt = 'Y \mum';% input('enter the xaxis label','s');
            yt = 'X \mum'; %input('enter the yaxis label','s');
            zt = 'Z \mum';
            label{1} = xlabel(xt);
            label{2} = ylabel(yt);
            label{3} = zlabel(zt);
            %set(gca,'YMinorTick','on','color',bcolor)
            ytickformat('%.1f')
            le{1} = 'plane 1'; %input('enter the legend','s');
            le{2} = 'plane 2';
            ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,1,[0 0])
            %errorbar(meanDisplacements(1,1:3),meanDisplacements(2,1:3),'.','color',[0 0 0],'MarkerSize',1)
            axis([10 30 25 60 0 ceil(max(r.Z))])
            legend off
            hold off
            view(0,90)
            savefile = 'vplanes4.tif';
            export_fig(vp,savefile,'-native');
        end
    end
end
