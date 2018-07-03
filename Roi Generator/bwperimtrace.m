function cellout = bwperimtrace(in,xlims,ylims)
% Traces the outer boundaries of regions of non-zero pixels in the input
% image (in), assuming a 4-connected neighbourhood
% ie. in [1 0 0; 0 1 0 ;0 0 0], the top-left and centre pixels will be
% considered separate objects. Output is a cell, with each element
% containing the coordinate list for the boundary of a separate object.
% Optional inputs xlims and ylims can be used to give the physical
% locations of the centre of the first and last pixels in the x and
% y directions respectively.

% Daniel Warren 10/11/2014
% Department of Oncology
% University of Oxford

if ~exist('xlims','var')
    xlims = 1:size(in,1);
end
if ~exist('ylims','var')
    ylims = 1:size(in,2);
end

if all(~in)
    cellout = cell(0);
    return;
end

in = in';

dx = (xlims(end)-xlims(1))/(size(in,1)-1);
dy = (ylims(end)-ylims(1))/(size(in,2)-1);

xpix = xlims(1):dx:xlims(2);
ypix = ylims(1):dy:ylims(2);

xbnd = [xpix-0.5*dx xpix(end)+0.5*dx];
ybnd = [ypix-0.5*dy ypix(end)+0.5*dy];

[indx, indy] = find(in);

points = zeros(8*length(indx),2);
%save('test.mat')
% Store line segments for each non-zero pixel's boundary (2 points for each
% segment, start and end), clockwise
x1 = xbnd(indx);
x2 = xbnd(indx+1);
y1 = ybnd(indy);
y2 = ybnd(indy+1);

points(1:8:end,1) = x1;
points(1:8:end,2) = y1;
points(2:8:end,1) = x2;
points(2:8:end,2) = y1;
points(3:8:end,:) = points(2:8:end,:);
points(4:8:end,1) = x2;
points(4:8:end,2) = y2;
points(5:8:end,:) = points(4:8:end,:);
points(6:8:end,1) = x1;
points(6:8:end,2) = y2;
points(7:8:end,:) = points(6:8:end,:);
points(8:8:end,:) = points(1:8:end,:);

% Remove opposing lines

midpts = 0.5*(points(1:2:end,:)+points(2:2:end,:));
dirs = (points(1:2:end,:)-points(2:2:end,:));

upinds = find(dirs(:,2) > 0);
downinds = find(dirs(:,2) < 0);
[~,upremv,downremv] = intersect(midpts(upinds,:),midpts(downinds,:),'rows');

rightinds = find(dirs(:,1) > 0);
leftinds = find(dirs(:,1) < 0);
[~,leftremv,rightremv] = intersect(midpts(leftinds,:),midpts(rightinds,:),'rows');

remove = [2*upinds(upremv);2*downinds(downremv);2*leftinds(leftremv);2*rightinds(rightremv)]-1;
remove = [remove;remove+1];

points = points(setdiff(1:end,remove),:);

% Reorder line segments, start new object if return to initial point
cellout = cell(0);
obj = 1;
finished = 0;
while ~finished
    i = 1;
    while i < size(points,1)
        [neighbind] = find(points(i+1,1) == points((i+2):end,1) & points(i+1,2) == points((i+2):end,2));
        neighbind(mod(neighbind,2) == 0) = [];
        if ~isempty(neighbind)
        neighbind = neighbind(1)+i+1;
        endind = (i+2):size(points,1);
        endind(endind == neighbind | endind == neighbind + 1) = [];
        points = [points(1:(i+1),:);points(neighbind+(0:1),:);points(endind,:)];
        end
        i = i+2;
        if all(points(i-1,:) == points(1,:))
            cellout{obj} = points(1:(i-1),:);
            if i > size(points,1)
                finished = 1;
            else
            points = points(i:end,:);
            obj = obj+1;
            end
            break
        end
    end
end

% Remove redundant points
for i = 1:length(cellout)
    cellout{i} = cellout{i}([1:2:(end-1) end],:);
    remove = [];
    for j = 2:(length(cellout{i})-1)
        prevdir = cellout{i}(j-1,:) - cellout{i}(j,:);
        nextdir = cellout{i}(j,:) - cellout{i}(j+1,:);
        if all(prevdir == nextdir)
            remove = [remove j];
        end
    end
    cellout{i} = cellout{i}(setdiff(1:end,remove),:);
end

end