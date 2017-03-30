clear all
close all
xyScale = 0.1625;
zScale = 0.4;
d = 3;
roiStack = getImages();
[nameTransFile,filePath] = uigetfile('*.tif','Select Transmitted Image for Overlay');
imageTrans = imread([filePath,nameTransFile]);
[nameBinary,filePath] = uigetfile('*.tif','Select Binary Image of Cell');
imageBinary = imread([filePath,nameBinary]);
roiStack=permute(roiStack, [2,1,3]);
imageSize(1,1) = size(roiStack,1)*xyScale;
imageSize(2,1) = size(roiStack,2)*xyScale;
imageSize(3,1) = size(roiStack,3)*zScale;

%%
sumImages = uint16(squeeze(max(permute(roiStack, [3,2,1]))));
sumImgScale = double(max(max(sumImages)))/(65536);
sumImages = uint16(sumImages/sumImgScale);
transImgScale = 65536/mean(prctile(imageTrans,95));
imageTrans = uint16(65536-double((imageTrans*transImgScale)));
transImgScale = 65536/mean(prctile(imageTrans,95));
imageTrans = uint16(double((imageTrans*(transImgScale/2))));
imshow(imageTrans)%invert (should make opaque objects brighter)
sumImages = sumImages+imageTrans; %combine dots and cells
sumImgScale = double(max(max(sumImages)))/65536;
sumImages = uint16(sumImages/sumImgScale);
imshow(sumImages,[]);
hold on
w = msgbox('Select a location with low displacements and double-click to continue');
                waitfor(w);
[~,sumBounds] = imcrop(sumImages);
close

sumBounds(1,3:4) = sumBounds(1,1:2) + sumBounds(1,3:4);
sumBounds(1,5:6) = (sumBounds(1,1:2) + sumBounds(1,3:4))/2;


x = size(roiStack,1);
y = size(roiStack,2);
z = size(roiStack,3);
%% Kilfoil Stack Filter
res=bpass3dMB(roiStack, [1 1 1], [12 12 12],[0 0]);

%% Kilfoil Object Detection 3D
masscut = mean(prctile(roiStack(:,:,size(roiStack,3)),95));
r=feature3dMB(res, d , [d d 10], [x y z],[1 1 1],d-1,masscut,.3); %
r(:,1:2) = r(:,1:2)*xyScale;
r(:,3) = r(:,3)*zScale;

%% Dots in Cell Region
rNDC = zeros(1,1);
for i = 1:size(r,1)
    %if it is under the cell
    if imageBinary(round(r(i,2)/xyScale),round(r(i,1)/xyScale))~=0     
        rNDC = cat(1,rNDC,i);
    end
end
rNDC(1,:) = [];

%% Dots in Cropped Region
sumBounds = sumBounds * xyScale;
sumIndX  = (r(:,1)>sumBounds(1,1) & r(:,1)<sumBounds(1,3));
sumIndY  = (r(:,2)>sumBounds(1,2) & r(:,2)<sumBounds(1,4));
sumIndXY = sumIndX .* sumIndY;
sumIndFinal = find(sumIndXY);
rND = (r(sumIndFinal,1:end));

% Find Z planes
h = histogram(rND(:,3),linspace(1,size(roiStack,3),size(roiStack,3)*3));
[~,zCenters] = find(imregionalmax(imgaussfilt(h.Values,3)));
zCenters = zCenters/3;
sumBounds(1,7) = zCenters(1,1);

% Dots in bottom Plane
for i =1:size(rND,1)
    rND(i,8) = find(min(abs(zCenters-rND(i,3)))==abs(zCenters-rND(i,3)));
end
bottomInd = rND(:,8)==1;
bottomIndFinal = find(bottomInd);
rNDB = rND(bottomIndFinal,1:end);

%Find Center dot (likely to have 4 equidistant neighbors)
for i = 1:size(rNDB,1)
differences(i,1) = rNDB(i,1)-sumBounds(1,5);
differences(i,2) = rNDB(i,2)-sumBounds(1,6);
differences(i,3) = rNDB(i,3)-sumBounds(1,7);
differences(i,4) = sqrt(differences(i,1)^2 + differences(i,2)^2 + differences(i,3)^2);
end
best = find(differences(:,4)==min(differences(:,4)));

k = 0;
count = 0;
while k == 0
%Find the 4 'equidistant' neighbors
clear differences
for i = 1:size(rNDB,1)
differences(i,1:3) = rNDB(i,1:3)-rNDB(best,1:3);
differences(i,4) = sqrt(differences(i,1)^2 + differences(i,2)^2 + differences(i,3)^2);
end
[sortedNew, sortedOrig] = sort(differences(:,4));
if std(sortedNew(2:5))<1
neighbors = rNDB(sortedOrig(2:5),:); 
neighbors(1:4,9) = sortedOrig(2:5);
k = 1;
elseif count > 20
    k = 1;
    disp('Could not find a suitable candidate for line fit')
else
newGuess = 2+round(3*rand());
best = sortedOrig(newGuess);
count = count +1;
end

end
close all
scatter3(neighbors(:,1),neighbors(:,2),neighbors(:,3))
hold on
scatter3(rNDB(best,1),rNDB(best,2),rNDB(best,3))
%%
% Pair off neighbors
clear differences
for i = 1:4
    for j = 1:4
    differences(j,1:3) = neighbors(i,1:3)-neighbors(j,1:3);    
    differences(j,4) = sqrt(differences(j,1)^2 + differences(j,2)^2 + differences(j,3)^2);
    end
    neighbors(i,10) = find(differences(:,4)==max(differences(:,4)));
end

used = zeros(1,4);
for i = 1:4
    if ismember(i,used)
    else
    dFit{i}(1,1:3) = neighbors(i,1:3);
    dFit{i}(2,1:3) = neighbors(neighbors(i,10),1:3);
    dFit{i}(3,1:3) = rNDB(best,1:3);
    used(i,1) = neighbors(i,10);
    end
end

%%

v1 = (dFit{1}(1,1:2) - dFit{1}(2,1:2))/2;
v2 = (dFit{2}(1,1:2) - dFit{2}(2,1:2))/2;

clear v1row
v1row = best;
for i = 1:2
    dv1 = rNDB(best,1:2);
while (dv1(1,1) < sumBounds(1,3) && dv1(1,1) > sumBounds(1,1) && dv1(1,2) < sumBounds(1,4) && dv1(1,2) > sumBounds(1,2)) == 1
   if i == 1
       clear differences
       dv1 = dv1 + v1;
       dv11(1:size(rNDB,1),1)=dv1(1,1);
       dv11(1:size(rNDB,1),2)=dv1(1,2);
       differences(:,1:2) = (rNDB(:,1:2) - dv11(:,1:2));
       differences(:,3) = sqrt(differences(:,1).^2+differences(:,2).^2);
       [dvSortNew,dvSortOrig] = sort(differences(:,3));
       if dvSortNew(1,1)<8
       v1row = cat(1,v1row,dvSortOrig(1,1));
       end
   else
       clear differences
       dv1 = dv1 + (v1*-1);
       dv11(1:size(rNDB,1),1)=dv1(1,1);
       dv11(1:size(rNDB,1),2)=dv1(1,2);
       differences(:,1:2) = (rNDB(:,1:2) - dv11(:,1:2));
       differences(:,3) = sqrt(differences(:,1).^2+differences(:,2).^2);
       [dvSortNew,dvSortOrig] = sort(differences(:,3));
       if dvSortNew(1,1)<8
       v1row = cat(1,v1row,dvSortOrig(1,1));
       end
   end
end
end

clear v2row
v2row = best;
for i = 1:2
    dv2 = rNDB(best,1:2);
while (dv2(1,1) < sumBounds(1,3) && dv2(1,1) > sumBounds(1,1) && dv2(1,2) < sumBounds(1,4) && dv2(1,2) > sumBounds(1,2)) == 1
   if i == 1
       clear differences
       dv2 = dv2 + v2;
       dv22(1:size(rNDB,1),1)=dv2(1,1);
       dv22(1:size(rNDB,1),2)=dv2(1,2);
       differences(:,1:2) = (rNDB(:,1:2) - dv22(:,1:2));
       differences(:,3) = sqrt(differences(:,1).^2+differences(:,2).^2);
       [dvSortNew,dvSortOrig] = sort(differences(:,3));
       if dvSortNew(1,1)<8
       v2row = cat(1,v2row,dvSortOrig(1,1));
       end
   else
       clear differences
       dv2 = dv2 + (v2*-1);
       dv22(1:size(rNDB,1),1)=dv2(1,1);
       dv22(1:size(rNDB,1),2)=dv2(1,2);
       differences(:,1:2) = (rNDB(:,1:2) - dv22(:,1:2));
       differences(:,3) = sqrt(differences(:,1).^2+differences(:,2).^2);
       [dvSortNew,dvSortOrig] = sort(differences(:,3));
       if dvSortNew(1,1)<8
       v2row = cat(1,v2row,dvSortOrig(1,1));    
       end
   end
end
end
v1row = unique(v1row);
v2row = unique(v2row);

[v1A,v1B] = fitLine3D(rNDB(v1row,1),rNDB(v1row,2),rNDB(v1row,3));
[v2A,v2B] = fitLine3D(rNDB(v2row,1),rNDB(v2row,2),rNDB(v2row,3));

for i=1:size(v1row,1)
    v1row(i,2) = norm(cross(v1B-v1A,rNDB(v1row(i,1),1:3)'-v1A))/norm(v1B-v1A);
end

for i=1:size(v2row,1)
     v2row(i,2) = norm(cross(v2B-v2A,rNDB(v2row(i,1),1:3)'-v2A))/norm(v2B-v2A);
end
v1mean = mean(v1row(:,2));
v2mean = mean(v2row(:,2));

if v2mean > v1mean
    rowV = (v1B-v1A)';
else
    rowV = (v2B-v2A)';
end



%%
[r,rows] = buildRows(r,rowV);

disp('done building rows')
%%
rowsNDC = rows;
for i = 1:size(rows,1)
    for j = 1:size(rows,2)
        if ismember(rows(i,j),rNDC) == 0
            rowsNDC(i,j) = 0;
        end
    end
    
    rowsNDC(i,find(rowsNDC(i,:)==0)) = max(rowsNDC(i,:));
    rowsNDCU(i,1:size(unique(rowsNDC(i,:)),2)) = unique(rowsNDC(i,:));
end

%%

clear rowsLines
for i = 1:size(rowsNDCU,1)
    n = nnz(rowsNDCU(i,:));
    if n >0
    [xyzFinal,rowP] = transLine3D(rowV,r(rowsNDCU(i,1:n),1:3),imageSize);
    rowsLines(1,1:3,i) = xyzFinal-rowV';
    rowsLines(2,1:3,i) = xyzFinal+rowV';
    end
end
disp('done fitting lines')
%%

figure
hold on
for i = 1:size(rows,1)
    n = nnz(rows(i,:));
    scatter3(r(rows(i,1:n),1),r(rows(i,1:n),2),r(rows(i,1:n),3))
end

for i = 1:size(rowsLines,3)
    plot3([rowsLines(1,1,i) rowsLines(2,1,i)],[rowsLines(1,2,i) rowsLines(2,2,i)],[rowsLines(1,3,i) rowsLines(2,3,i)])
end

scatter3(0,0,0)
%%
rRef = zeros(size(r,1),3);
for i = 1:size(rowsLines,3)
    if nnz(rowsNDCU(i,:))>0
    n = nnz(rows(i,:));
    [xy,distance,t_a] = distance2curve(rowsLines(:,:,i),r(rows(i,1:n),1:3));
    for j = 1:size(xy,1)
        rRef(rows(i,j),1:3) = xy(j,1:3);
    end
    else
        for j = 1:nnz(rows(i,:))
        rRef(rows(i,j),1:3) = r(rows(i,j),1:3);
        end
    end
end
rDisp = r(:,1:3)-rRef(:,1:3);
%%
rDispFilt = rDisp;
rDispFilt((rDispFilt<.4 & rDispFilt>-.4))=0;

meshRes = 2;
%[xx,yy,zz] = meshgrid(1:meshRes:size(roiStack,1),1:meshRes:size(roiStack,2),1:1:size(roiStack,3));
[xx,yy,zz] = meshgrid(1:meshRes:imageSize(1,1),1:meshRes:imageSize(2,1),min(r(:,3)):0.2:max(r(:,3)));
xq = double([xx(:) yy(:) zz(:)]);
vq = griddatan(rRef(:,1:3),rDispFilt(:,3),xq);
vq = reshape(vq,size(xx));
vq(isnan(vq)) = min(min(min(vq)));

vq2 = vq;
vq2 = vq2+abs(min(min(min(vq2))));
vq2Scale = double(65000/max(max(max(vq2))));
vq2 = uint16(vq2Scale*vq2);
ShowStack(vq2,1,1)


%%
% close all
% scatter3(rNDB(v1row(:,1),1),rNDB(v1row(:,1),2),rNDB(v1row(:,1),3))
% hold on
% scatter3(rNDB(v2row(:,1),1),rNDB(v2row(:,1),2),rNDB(v2row(:,1),3))
% scatter3(rNDB(best,1),rNDB(best,2),rNDB(best,3))
% scatter3(neighbors(:,1),neighbors(:,2),neighbors(:,3))
% plot3([v1A(1,1),v1B(1,1)],[v1A(2,1),v1B(2,1)],[v1A(3,1),v1B(3,1)])
% plot3([v2A(1,1),v2B(1,1)],[v2A(2,1),v2B(2,1)],[v2A(3,1),v2B(3,1)])

%%
% %%
% clear mask
% mask = find(r(:,3)<20);
% layer = fit([r(mask,1),r(mask,2)],r(mask,3),'lowess','Span',0.005);
% %%
% close all
% %scatter3(r(mask,1),r(mask,2),r(mask,3));
% 
% plot(layer)
% hold on
% scatter3(0,0,0);