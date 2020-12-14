function [plane,planesLocFiltList] = placePlanes(r,plane,Surface2)
%%    
if isempty(Surface2)~=1
    for j = 1:size(plane.final,2)
        for i = 1:nnz(plane.final(:,j))
            plane.pos(i,j) = (feval(Surface2,r.X(plane.final(i,j)),r.Y(plane.final(i,j)))) - r.Z(plane.final(i,j));
        end
    end
    
else
    for j = 1:size(plane.final,2)
        for i = 1:nnz(plane.final(:,j))
            plane.pos(i,j) = max(r.Z) - r.Z(plane.final(i,j));
        end
    end    
end


plane.pos(plane.pos==0)=nan;
plane.loc2 = mean(plane.pos,'omitnan');
if min(plane.loc2) < 0
    plane.loc2 = plane.loc2 + abs(min(plane.loc2));
end

planesLocFilt = find(plane.loc2>4.5);
planesLocFiltList = plane.final(:,planesLocFilt(1,1));
if size(planesLocFilt,2)>1
    for i = 2:size(planesLocFilt,2)
        planesLocFiltList = cat(1,planesLocFiltList,plane.final(:,planesLocFilt(1,i)));
    end
end
planesLocFiltList(planesLocFiltList==0) = [];



%Old Approach
% fS = open('Shear Mat Files\Surface.mat');
% try
%     fitSurface = fS.fitSurface{3};
% catch
%     xtemp = [1;2;3;4;5];
%     ytemp = [5;8;3;5;7];
%     ztemp = size(res,3)*raw.dataKey(10,1)*ones(5,1);
%     fitSurface = fit([xtemp,ytemp],ztemp,'lowess','Span',0.1);
% end