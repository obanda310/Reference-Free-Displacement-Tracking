function [r,rows,plane] =VarUpkeep(r,rows,plane)
%This function will use the current instance of the the RawData3D class
%object to update related classes' information. This should always solve
%any upkeep related issues as long as the RawData3D is treated as the
%'master' data set.

% Before beginning, determine what stage of the program we are at.
r.State = 0;
if size(r.rVa,1) > 0
    r.State = 2;
elseif size(r.rX,1)>0
    r.State = 1;
end


%% First Check for obvious errors and correct
%Error type 1: 2 objects in same column have same plane

%Solution: determine which object better matches the rest of the
%column(if data is available), otherwise arbitrarily assign a new column to
%one of the objects

for i = 1:max(r.col)
    idx = r.col == i;
    if nnz(r.plane(idx)) > nnz(unique(r.plane(idx)))
        tplane = mode(r.plane(idx));
        probIdx = find(r.col ==i & r.plane==tplane);
        goodIdx = find(r.col ==i & r.plane~=tplane);
        if size(goodIdx,1) > 0 %if other column data exists
            %tdist = pdist(r.r(probIdx,1:3)); %may add additional logic
            %using this variable later
            keep = dsearchn(r.r(probIdx,1:2),r.r(goodIdx(1,1),1:2));
            probIdx(keep,:) = [];
            for j = 1:size(probIdx,1)
                r.col(probIdx(j)) = max(r.col)+1;
            end
        else % otherwise, arbitrarily reassign
            probIdx(1,:) = [];
            for j = 1:size(probIdx,1)
                r.col(probIdx(j)) = max(r.col)+1;
            end
        end
    end
end

%Error type 2: Multiple detections in close proximity
%Solution: If the objects are not in the top plane merge them and delete
%one. Since choosing which object to keep will require objects to be sorted
%already, only perform this step if rState is greater than 0.
removeIdx = [];
if r.State>0
    for i = 2:max(r.planeGroup)
        tidx = find(r.planeGroup == i);
        pPts = rangesearch(r.r(tidx,1:2),r.r(tidx,1:2),.3);
        for j = 1:size(pPts,1)
            if nnz(pPts{j})>1
                idx = sort(pPts{j}(:));
                picked = 0;
                for k = 1:nnz(idx)
                    if r.colUp(tidx(idx(k))) == 0 && picked == 0
                        removeIdx = cat(1,tidx(idx(k)),removeIdx);
                        picked = 1;
                    end
                end
            end
        end
    end
end
removeIdx = unique(removeIdx);

r = removeIndex(r,removeIdx);

%Error type 3: no row or col assignment and XY position is near image
%border
%Solution: Delete it.
borderwidth = 6;
posIdx = r.r(:,1)<=borderwidth | r.r(:,2)<=borderwidth | r.s(2,2)-r.r(:,1)<=borderwidth | r.s(1,2)-r.r(:,2)<=borderwidth;
rcCheckIdx = (r.col==0 & r.row==0);
removeIdx = posIdx & rcCheckIdx;

r = removeIndex(r,removeIdx);

%Error type 4: Column assignment exists and is in a VPlane, but no row exists.
%Solution: Assign row as the most common non-zero row in that plane within
%the vplane.
if r.State>0
    tidx = find(r.col>0 & r.vplane>0 & r.row==0);
    for i = 1:nnz(tidx)
        r.row(tidx(i)) = mode(r.row((r.vplane==r.vplane(tidx(i)))&(r.plane==r.plane(i))&(r.row>0)));
    end
end
%% Sort the Data

%First, sort all the data in obj. Start at the origin, and go
%row by row (using row vector), keep columns together.
if r.State == 0
    if abs(rows.V(1,1))<abs(rows.V(2,1))
        [~,order] = sortrows(r.r,[1,2]);
    else
        [~,order] = sortrows(r.r,[2,1]);
    end
elseif r.State == 1 %If reference info exists, sort by that instead
    if abs(rows.V(1,1))<abs(rows.V(2,1))
        [~,order] = sortrows([r.rX,r.rY,r.vplane],[3,1,2]);
    else
        [~,order] = sortrows([r.rX,r.rY,r.vplane],[3,2,1]);
    end
elseif r.State == 2 %If reference info exists, sort by that instead (currently redundant code, may add steps later)
    if abs(rows.V(1,1))<abs(rows.V(2,1))
        [~,order] = sortrows([r.rX,r.rY,r.vplane],[3,1,2]);
    else
        [~,order] = sortrows([r.rX,r.rY,r.vplane],[3,2,1]);
    end
end

clear tidx
tidx(:,1) = (1:1:size(r.X,1));
temp(:,1) = r.col(order);
temp(:,2) = r.col(order);
temp(:,3) = tidx(order);

cc = 1; %current column count
for i = 1:size(temp,1)
    if temp(i,1) > cc
        temp(temp(:,1)==cc,1) = max(temp(:,1))+1; % make cc idx available for safe assignment
        temp(temp(:,1)==temp(i,1)) = cc; % assign all objects with current column idx to cc idx
        cc = cc+1;
    end
end

temp = sortrows(temp,3);

r.r(:,8) = temp(:,1);


[r.r,orderF] = sortrows(r.r,[8,3]);
%% Sort all other elements of class object by new order

r.X = (r.X(orderF));
r.Y = (r.Y(orderF));
r.Z = (r.Z(orderF));
r.Mass = (r.Mass(orderF));
r.MaxI = (r.MaxI(orderF));
r.rg = (r.rg(orderF));
r.PctAbove = (r.PctAbove(orderF));
r.row = (r.row(orderF));
r.col = r.r(:,8);
r.colSize = (r.colSize(orderF));
r.plane = (r.plane(orderF));
r.planeGroup = (r.planeGroup(orderF));
r.Dl = (r.Dl(orderF));
r.NDl = (r.NDl(orderF));

if r.State>0 % these only work after vplanes have been created
    r.rX = (r.rX(orderF));
    r.rY = (r.rY(orderF));
    r.rZ = (r.rZ(orderF));
    r.colUp = (r.colUp(orderF));
    r.colDown = (r.colDown(orderF));
    r.rowLeft = (r.rowLeft(orderF));
    r.rowRight = (r.rowRight(orderF));
    r.colCheck = (r.colCheck(orderF));
    r.vplane = (r.vplane(orderF));
end

if r.State>1 % this only works after row fits have been created
    r.rVa = (r.rVa(orderF,:));
    r.rVb = (r.rVb(orderF,:));
    r.dS = (r.dS(orderF,:));
    r.dX = (r.dX(orderF,:));
    r.dY = (r.dY(orderF,:));
    r.dZ = (r.dZ(orderF,:));
    
end
%
r.D = find(r.Dl);
r.ND = find(r.NDl);

%% Rebuild Rows class object
%first remove empty row indices
for i = 1:max(r.row)
    if nnz(r.row==i) ==0
        r.row(r.row>i) = r.row(r.row>i)-1;
    end
end

rows.m = zeros(max(r.row),1);
rows.p = zeros(max(r.row),1);
rows.NDm = zeros(max(r.row),1);
rows.s = zeros(max(r.row),1);
for i = 1:size(r.row,1)
    if r.row(i)>0
        if rows.m(r.row(i),1) == 0
            rows.m(r.row(i),1) = i;
        else
            rows.m(r.row(i),nnz(rows.m(r.row(i),:))+1) = i;
        end
        
        if r.NDl(i) == 1 && rows.NDm(r.row(i),1) == 0
            rows.NDm(r.row(i),1) = i;
        elseif r.NDl(i) == 1
            rows.NDm(r.row(i),nnz(rows.NDm(r.row(i),:))+1) = i;
        end
    end
end
for i = 1:size(rows.m,1)
    rows.p(i) = median(r.plane(rows.m(i,1:nnz(rows.m(i,:)))));
    rows.s(i) = nnz(rows.m(i,:));
end

%% Rebuild Planes Class Object
plane.final = zeros(size(plane.final));
for i = 1:size(r.row,1)
    if r.plane(i)>0
        if plane.final(1,r.plane(i)) == 0
            plane.final(1,r.plane(i)) = i;
        else
            plane.final(nnz(plane.final(:,r.plane(i)))+1,r.plane(i)) = i;
        end
    end
end

    function r = removeIndex(r,removeIdx)
        r.r(removeIdx,:) = [];
        r.X(removeIdx,:) = [];
        r.Y(removeIdx,:) = [];
        r.Z(removeIdx,:) = [];
        r.Mass(removeIdx,:) = [];
        r.MaxI(removeIdx,:) = [];
        r.rg(removeIdx,:) = [];
        r.PctAbove(removeIdx,:) = [];
        r.row(removeIdx,:) = [];
        r.col(removeIdx,:) = [];
        r.colSize(removeIdx,:) = [];
        r.plane(removeIdx,:) = [];
        r.planeGroup(removeIdx,:) = [];
        r.Dl(removeIdx,:) = [];
        r.Dl = r.Dl>0; %convert to logical
        r.NDl(removeIdx,:) = [];
        r.NDl = r.NDl>0; %convert to logical
        
        if r.State>0 % these only work after vplanes have been created
            r.rX(removeIdx,:) = [];
            r.rY(removeIdx,:) = [];
            r.rZ(removeIdx,:) = [];
            r.colUp(removeIdx,:) = [];
            r.colDown(removeIdx,:) = [];
            r.rowLeft(removeIdx,:) = [];
            r.rowRight(removeIdx,:) = [];
            r.vplane(removeIdx,:) = [];
            r.colCheck(removeIdx,:) = [];
        end
        
        if r.State>1 % this only works after row fits have been created
            r.rVa(removeIdx,:) = [];
            r.rVb(removeIdx,:) = [];
            r.dS(removeIdx,:) = [];
            r.dX(removeIdx,:) = [];
            r.dY(removeIdx,:) = [];
            r.dZ(removeIdx,:) = [];
        end
        
        r.l = size(r.X,1);
    end
end