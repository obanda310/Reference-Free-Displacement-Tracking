function [r,VPlanes,rows,plane] = RefreshVars(r,rows,plane)
%% Clean all major class variables
[r,rows,plane] = VarUpkeep(r,rows,plane);
%disp(['Refreshing Variables'])
[VPlanes,r] = VertPlaneData(r,rows);
VPlanes = removeDuplicates(VPlanes);
r = updateVPlane(r,VPlanes,0);
[VPlanes,r] = ProfileColumns(VPlanes,r,plane);
[VPlanes,r] = BuildGrids(VPlanes,r,rows);
if r.State>1
    VPlanes = BuildDispGrids(VPlanes,r); %Added step now that displacements are known
end
r = updateVPlane(r,VPlanes,0);
[VPlanes,r] = FillHoles(VPlanes,r,1);
[rows,r]=formatRows(rows,plane,r);

try
    %disp(['done refreshing at ' num2str(toc) ' seconds'])
catch
end
end