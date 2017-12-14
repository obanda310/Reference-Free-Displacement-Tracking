function [topSurface ,fitSurface] = findSurface(shear,cm2,imageArea,imageBorders,dataKey)

%clear topSurface topSurfaceCand bLimits fitSurface

%find all pillars with max disp below threshold
topSurfaceCand = cm2(1:nnz(cm2(:,1)),1);
for i = 2:(shear.cmCutoff)
    topSurfaceCand = cat(1,topSurfaceCand,cm2(1:nnz(cm2(:,i)),i));
end


%get rid of pillars beneath the cell (might have normal displacement)
topSurface = zeros(1,1);

for i = 1:size(topSurfaceCand,1)
    %if it is under the cell
    if imageArea(round(shear.rawY1(topSurfaceCand(i,1))/dataKey(9,1)),round(shear.rawX1(topSurfaceCand(i,1))/dataKey(9,1)))~=0 && imageArea(round(shear.lastY(topSurfaceCand(i,1))/dataKey(9,1)),round(shear.lastX(topSurfaceCand(i,1))/dataKey(9,1)))~=0
        %if it is near the image border
        if imageBorders(round(shear.rawY1(topSurfaceCand(i,1))/dataKey(9,1)),round(shear.rawX1(topSurfaceCand(i,1))/dataKey(9,1)))~=0 && imageBorders(round(shear.lastY(topSurfaceCand(i,1))/dataKey(9,1)),round(shear.lastX(topSurfaceCand(i,1))/dataKey(9,1)))~=0
            topSurface = cat(1,topSurface,topSurfaceCand(i,1));
        end
    end
end
%shift cells up 1 to get rid of initial zero
topSurface(1,:) = [];

for i = 1:size(topSurface,1)
    topSurface(i,2) = shear.rawX(shear.lastFrame(topSurface(i,1)),topSurface(i,1));
    topSurface(i,3) = shear.rawY(shear.lastFrame(topSurface(i,1)),topSurface(i,1));
    topSurface(i,4) = shear.lastFrame(topSurface(i,1));
end

%remove outliers (typically incomplete pillars at the edges)
ub = mean(topSurface(:,4)) + 2*std(topSurface(:,4));
lb = mean(topSurface(:,4)) - 2*std(topSurface(:,4));
topSurface([find(topSurface(:,4)>ub | topSurface(:,4)<lb)],:) = [];


% fitSurface{1} = fit([topSurface(:,3),topSurface(:,2)],topSurface(:,4),'poly11');
% fitSurface{2} = fit([topSurface(:,3),topSurface(:,2)],topSurface(:,4),'lowess','Span',0.1);
fitSurface{3} = fit([topSurface(:,3),topSurface(:,2)],topSurface(:,4)*dataKey(10,1),'lowess','Span',0.1);