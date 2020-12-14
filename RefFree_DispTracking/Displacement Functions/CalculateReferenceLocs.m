function [r,m3,image] = CalculateReferenceLocs(r,rows,plane,VPlanes,image,raw,planesLocFiltList)
disp('Updating Reference Locations')
% 1st METHOD FOR MEASURING DISPLACEMENTS
%disp('Fitting Rows with Independent Slopes - Method 1')
m1 = DispData3D;
m1 = method1fit(m1,r,rows);
%disp(['done Fitting Rows with Independent Slopes - Method 1 at ' num2str(toc) ' seconds'])

% 2nd METHOD FOR MEASURING DISPLACEMENTS V2
%disp('Fitting Rows with Row Slope - Method 2')
m2 = DispData3D;
m2 = method2fit(m2,m1,r,raw,image.ADil,rows);
%disp(['done Fitting Rows with Row Slope - Method 2 at ' num2str(toc) ' seconds'])

% 3RD METHOD FOR FITTING ROWS (METHOD 1 AND 2 COMBINED)
%disp('Picking Best Case Row Fit - Method 3')
m3 = DispData3D;
m3 = method3fit(m3,m1,m2,rows,raw);
% Calculate Displacements from Fit Lines Method 3
m3 = calcDisp(m3,r,VPlanes); %based on fits of non-deformed markers
m3 = calcDispSC(m3,r,VPlanes); %includes shift correction in reference approximation. ****Not in use currently****
m3.disp = r.r(:,1:3)-m3.ref(:,1:3);
% Calculate Average Displacement per Row per Plane Method 3
m3 = dispStats(m3,plane,rows,r,planesLocFiltList);
% Filter out noise in Displacements Method 1
m3 = dispNoise(m3,r,planesLocFiltList,image,raw.dataKey(9,1));
image.imgNBds = imageNoiseBounds(r.r,m3,image,raw.dataKey(9,1));
%disp(['done Picking Best Case Row Fit - Method 3 at ' num2str(toc) ' seconds'])
r = updateReference(r,m3);
disp(['done Updating Reference Locations at ' num2str(toc) ' seconds'])