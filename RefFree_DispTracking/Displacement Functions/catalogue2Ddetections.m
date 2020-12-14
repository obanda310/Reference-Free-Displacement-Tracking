function [xyz2Dfinal,xyz2D,lub] = catalogue2Ddetections(image,raw,r)

% Sub-pixel Object Detection Using Kilfoil 'feature2D' Function
clear xyz2D
%[r] = feature2D(img,lambda,w,masscut,Imin)
for i = 1:size(image.MaskStack,3)
    clear temp
    currentImg = image.MaskStack(:,:,i);
    temp = feature2D(currentImg,1,2,2000,50);
    if i == 1
        temp(:,6) = i;
        xyz2D = temp;
    else
        temp(:,6) = i;
        xyz2D = cat(1,xyz2D,temp);
    end
end

% Here we will link objects in the Z direction
disp('Linking 3D Detections Between Planes.')
% Set a maximum linking distance in microns that any object can still be
% considered part of a pillar. Smaller values will speed up code.
maxLinkDistance = 1.2; % !!**value should be smaller than array spacing**!!
maxLD = maxLinkDistance;% /raw.dataKey(9,1); <--if using pixels, use this conversion
disp(['Maximum allowable XY Translation (Microns): ',num2str(maxLinkDistance)])
% Set a maximum number of frames to look for a linked object before giving
% up (maxJumpDistance)
maxJD = 1;
% Linking using trackmem function (kilfoil code)
disp('Linking Pillars')
[lub] = trackmem(xyz2D,maxLD,2,0,maxJD);

xyz2Ds = xyz2D;
xyz2Ds(:,3:5) = [];
xyz2Ds(:,1) = xyz2D(:,2);
xyz2Ds(:,2) = xyz2D(:,1);
xyz2Ds(:,1:2) = xyz2Ds(:,1:2)*raw.dataKey(9,1); %convert to microns
xyz2Ds(:,3) = xyz2Ds(:,3)*raw.dataKey(10,1); %convert to microns
%% Object Detection Plot
detections = figure;
imshow(image.Black)
hold on
map = brewermap(max(xyz2D(:,6)),'*Spectral');

for i = 1:max(xyz2D(:,6)) %size(subpixMaxima,3)-10
    [tempFrame,~] = find(xyz2D(:,6) == i);
    scatter3(xyz2D(tempFrame,1),xyz2D(tempFrame,2),xyz2D(tempFrame,6),'.','SizeData',100,'MarkerFaceColor',map(i,1:3),'MarkerEdgeColor',map(i,1:3))
end
hold off
filePath = cd;
savefile = [filePath '\Tracking_Unlinked Detections.tif'];
export_fig(detections,savefile,'-native');


%% Link 3D detections to 2D Detections and Apply Corrections
% 2D Centroids (in XY) can provide locations of 'incomplete' features near the
% surface, where deformations are highest, and where 3D object detection can fail
% to detect objects. Unfortunately, 2D info is subject to additional flaws,
% such as tilted ellipsoids (caused by imperfect laser alignment), or
% tilted columns (caused by misaligned patterning plane, imaging plane,
% and hydrogel surface. This section will attempt to use 3D info as a
% guide to catalogue 2D info.

% Preserve original 2D list
xyz2Di = xyz2Ds;
%%
xyz2Ds = xyz2Di; %use to reset
% Creat a final list with 3D assignments
xyz2Dfinal = xyz2Ds;

% First link 3D objects to the 2D objects within a maximum allowable
% distance
maxLD = 1.5; %should be less than half of the Z spacing during patterning
clear xyz2Dci xyz2Ddi rDe
rDe = delaunayn(double([r.X, r.Y, r.Z]));
[xyz2Dci,xyz2Ddi] = dsearchn(double([r.X, r.Y, r.Z]),rDe,xyz2Ds); %indices and distance of closest 2D detections
xyz2Dci(xyz2Ddi(:,1)>maxLD,1) = 0;
xyz2Ddi(xyz2Ddi(:,1)>maxLD,1) = 0;
xyz2Dfinal(:,4) = xyz2Dci(:,1);
xyz2Dfinal(:,5) = xyz2Ddi(:,1);
xyz2Dfinal(xyz2Dfinal(:,4)>0,6) = r.Z(xyz2Dfinal(xyz2Dfinal(:,4)>0,4))-xyz2Dfinal(xyz2Dfinal(:,4)>0,3);
xyz2Dfinal = sortrows(xyz2Dfinal,[4 6]);
toc

lowCO = raw.dataKey(10,1) - maxLD;
highCO = maxLD - raw.dataKey(10,1);
maxLD = 1;
go =1;
while go == 1
    clear pEnds
    sizeBef = size(xyz2Ds,1);
    xyz2Ds = xyz2Dfinal(xyz2Dfinal(:,4)==0,1:3);
    sizeAft = size(xyz2Ds,1);
    if sizeBef == sizeAft
        break
    end
    pEnds = cat(1,xyz2Dfinal(xyz2Dfinal(:,6)<lowCO,1:3),xyz2Dfinal(xyz2Dfinal(:,6)>highCO,1:3));
    pIdx = cat(1,xyz2Dfinal(xyz2Dfinal(:,6)<lowCO,4),xyz2Dfinal(xyz2Dfinal(:,6)>highCO,4));
    try
        rDe = delaunayn(pEnds);
    catch
        break
    end
    [xyz2Dci,xyz2Ddi] = dsearchn(pEnds,rDe,xyz2Ds);
    xyz2Dci(xyz2Ddi(:,1)>maxLD,1) = 0;
    xyz2Ddi(xyz2Ddi(:,1)>maxLD,1) = 0;
    xyz2Didx = zeros(size(xyz2Dci));
    xyz2Didx(xyz2Dci(:,1)>0,1) = pIdx(xyz2Dci(xyz2Dci(:,1)>0),1);
    xyz2Dfinal(1:size(xyz2Dci,1),4) = xyz2Didx(:,1);
    xyz2Dfinal(1:size(xyz2Ddi,1),5) = xyz2Ddi(:,1);
    xyz2Dfinal(xyz2Dfinal(:,4)>0,6) = r.Z(xyz2Dfinal(xyz2Dfinal(:,4)>0,4),1)-xyz2Dfinal(xyz2Dfinal(:,4)>0,3);
    xyz2Dfinal = sortrows(xyz2Dfinal,[4 6]);
    lowCO = lowCO - raw.dataKey(10,1);
    highCO = highCO + raw.dataKey(10,1);
    toc
end
disp('It ended')
%%
figure
scatter3(xyz2Dfinal(xyz2Dfinal(:,4)==0,1),xyz2Dfinal(xyz2Dfinal(:,4)==0,2),xyz2Dfinal(xyz2Dfinal(:,4)==0,3))
% Then use the farthest assigned 2D detections to continue assignments
% until no more assignments are possible

% Iteratively assign a 3D object index to a 2D object index if it is
% within the allowable distance, and remove that 2D object from the search
% list. Then repeat the search
%% Calculate average slope of pillars in Z

for i = 1:size(r.X,1)
    if r.colSize(i) >1 && ismember(i,r.ND)
        clear tempcol
        tempcol(:,1:3) = r.r(r.col(:)==r.col(i),1:3);
        sortrows(tempcol,3);
        Xs(i) = tempcol(1,1)-tempcol(end,1);
        Ys(i) = tempcol(1,2)-tempcol(end,2);
        Zs(i) = tempcol(1,3)-tempcol(end,3);
    end
end
colV(1) = mean(Xs);
colV(2) = mean(Ys);
colV(3) = mean(Zs);
colV = colV/min(colV);
