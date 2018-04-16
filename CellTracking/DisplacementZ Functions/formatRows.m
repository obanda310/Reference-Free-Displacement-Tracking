function [rows,rowPlanes,rowPlanesIdx,rowsNDCU,r]=formatRows(rows,plane,r)

% Separate Rows by Plane
clear rowPlanes
for i = 1:size(plane.final,2)
    for j = 1:size(rows,1)
        rowPlanes(j,1:size(intersect(rows(j,:),plane.final(1:nnz(plane.final(:,i)),i)),1),i) = intersect(rows(j,:),plane.final(1:nnz(plane.final(:,i)),i));
    end
end
%%
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
rows = newRows;

% Label Objects in 'r' with Their Row Number
for i = 1:size(rows,1)
    for j = 1:size(rows,2)
        if rows(i,j)>0
            r.row(rows(i,j)) = i;
        end
    end
end

% Determine non-deformed members of each row
clear rowsNDCU rowsNDC
rowsNDC = rows;
for i = 1:size(rows,1)
    for j = 1:size(rows,2)
        if ismember(rows(i,j),r.ND) == 0
            rowsNDC(i,j) = 0;
        end
    end
    rowsNDC(i,find(rowsNDC(i,:)==0)) = max(rowsNDC(i,:));
    rowsNDCU(i,1:size(unique(rowsNDC(i,:)),2)) = unique(rowsNDC(i,:));
end

end