function [rDisp1Mean,rDisp1MeanPlanes,rDisp1StdPlanes,rDisp1MeanTotal,rDisp1StdTotal] = quickRowStats(rDisp1,planesFinal,rowPlanesIdx,rNDC,rowsNDCU,planesLocFiltList)
%% Calculate Average Displacement per Row per Plane Method 2
for i = 1:size(rowPlanesIdx,1)
    for j = 1:(rowPlanesIdx(i,2)-rowPlanesIdx(i,1))+1
        rDisp1Mean(j,i) = mean(abs(rDisp1(rowsNDCU(j+(rowPlanesIdx(i,1)-1),rowsNDCU(j+(rowPlanesIdx(i,1)-1),:)>0),3)));
    end
    rDisp1MeanPlanes(1,i) = mean(abs(rDisp1(intersect(rNDC,planesFinal(1:nnz(planesFinal(:,i)),i)),3)));
    rDisp1StdPlanes(1,i) = std(abs(rDisp1(intersect(rNDC,planesFinal(1:nnz(planesFinal(:,i)),i)),3)));
end
rDisp1MeanTotal = mean(abs(rDisp1(intersect(planesLocFiltList(:,1),rNDC),3)));
rDisp1StdTotal = std(abs(rDisp1(intersect(planesLocFiltList(:,1),rNDC),3)));
end
