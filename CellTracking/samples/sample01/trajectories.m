close all; clear; clc;

% for u(ser)
[name,path] = uigetfile('*.xlsx');
file = [path,name];

% % for me
% file = 'F:\trajectories.xlsx';
[num,txt,raw] = xlsread(file);

traj = num(:,2);
numTraj = max(traj); % number of objects

for i = 1:numTraj
    obj = num(traj==i,:);
    numFrames = size(obj,1);
    for j = 1:numFrames
        x_0(i,j) = obj(1,4);
        y_0(i,j) = obj(1,5);
        x_diff(i,j) = obj(j,4) - obj(1,4);
        y_diff(i,j) = obj(j,5) - obj(1,5);
    end
end
%%
% for k = 1:23        % number of z-slices
%     scale = 2.5;
%     data = figure('units','normalized','outerposition',[0 0 1 1]);
%     h = quiver(x_0(:,k),512-y_0(:,k),x_diff(:,k),-y_diff(:,k),0);
%     hU = get(h,'UData');
%     hV = get(h,'VData');
%     set(h,'UData',scale*hU,'VData',scale*hV)
%     axis(gca,'tight')
%     savefile = sprintf('Trajectory between 1st Frame and Frame %u.tif',k);
%     print(savefile,'-dtiff','-r1800')
%     close
% end

%%
objAll = [x_0(:,1),512-y_0(:,1)];
rows = zeros(2,2);
numObj = length(x_0);
row = 1;
col = 1;
% distances = cell(numObj,3);
for i = 1:numObj
    thisObj = [x_0(i,1),512-y_0(i,1)];
    objCount = (1:numObj)';
    distances = sqrt(sum(bsxfun(@minus,objAll,thisObj).^2,2));
    objData = [objCount,distances];
    sortDist = sortrows(objData,2);
    neighbors(1,1:4) = sortDist(2:5,1);
    nDistance = sortDist(2:5,2);
    for j = 1:4
        n(j,1:2) = objAll(neighbors(1,j),:);
        n(j,3) = neighbors(j);
    end
    a = bsxfun(@minus,n(:,1:2),thisObj);
    b = mean([min(abs(a)),max(abs(a))]);
    goodNeighbors = n(abs(a(:,2)) < b & a(:,1) > 0,:);
    if size(goodNeighbors,1) > 1
        sortN = sortrows(goodNeighbors,1);
        rightNeighbor = sortN(1,:);
    else
        rightNeighbor = goodNeighbors;
    end
    
    if isempty(rightNeighbor)
        row = row + 1;
        col = 1;
    else
        if isempty(rightNeighbor(3)) == 0
            [trow,tcol] = find(rows(row,col),rightNeighbor(3));
            if isempty(trow) ~= 0
                rows(trow,end+1) = i;
            else
                rows(row,end+1) = i;
            end
        end
    end
end

%%
plot(objAll(:,1),objAll(:,2),'o')
hold on
plot(objAll(418,1),objAll(418,2),'ro')