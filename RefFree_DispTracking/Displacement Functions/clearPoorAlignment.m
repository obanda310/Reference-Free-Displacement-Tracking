function [r,VPlanes,rows,plane] = clearPoorAlignment(r,VPlanes,rows,plane)

%% Find columns where successive vector orientation doesn't match.
% (should be unidirectional and with increasing magnitude.)

% Find large deformations and verify lower objects in column match trajectory
idx = find(r.dS>.4 & r.rX~=0);
alignment = zeros(size(idx));
for i = 1:nnz(idx)
    tvec1 = [r.dX(idx(i)),r.dY(idx(i))]./norm([r.dX(idx(i)),r.dY(idx(i))]);
    if r.colDown(idx(i))>0
        tvec2 = [r.dX(r.colDown(idx(i))),r.dY(r.colDown(idx(i)))]./norm([r.dX(r.colDown(idx(i))),r.dY(r.colDown(idx(i)))]);
        if r.dS(r.colDown(idx(i))) > .2
            alignment(i,1) = dot((tvec1),(tvec2)./norm(tvec2));
        else
            alignment(i,1) = NaN;
        end
    end
    if r.rowRight(idx(i))>0
        tvec2 = [r.dX(r.rowRight(idx(i))),r.dY(r.rowRight(idx(i)))]./norm([r.dX(r.rowRight(idx(i))),r.dY(r.rowRight(idx(i)))]);        
        if r.dS(r.rowRight(idx(i))) > .5
            alignment(i,2) = dot((tvec1),(tvec2)./norm(tvec2));
        else
            alignment(i,2) = NaN;
        end
    end
    if r.rowLeft(idx(i))>0
        tvec2 = [r.dX(r.rowLeft(idx(i))),r.dY(r.rowLeft(idx(i)))]./norm([r.dX(r.rowLeft(idx(i))),r.dY(r.rowLeft(idx(i)))]);       
        if r.dS(r.rowLeft(idx(i))) > .5
            alignment(i,3) = dot((tvec1),(tvec2)./norm(tvec2));
        else
            alignment(i,3) = NaN;
        end
    end
    
end
alignment(alignment==0)=NaN;
tidx2 = idx(min(alignment,[],2)<0.75);
r = resetData(r,tidx2);
VPlanes = resetGridObject(VPlanes,tidx2);
% Clean all major class variables
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);