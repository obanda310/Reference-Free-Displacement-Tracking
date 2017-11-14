%function dispZfunc(directory)
% if nargin ==1
% cd(directory);
% end
clear all
close all
%%
set(0,'defaultfigurecolor',[1 1 1])
tic

%%
files = dir('*.tif'); %Check Directory for default filenames
matFiles = dir('*.mat');



xyScale = 0.1625;
zScale = 0.4;
d = 3;

if size(matFiles,1)>=1
    for k = 1:length(matFiles)
        current=matFiles(k).name;
        if size(current,2)>12
            check(k)=strcmp(current(end-12:end),'StackName.mat');
        else
            check(k) = 0;
        end
    end
end
try
    if size(find(check),2) >0
        load('StackName.mat','StackName')
        roiStack = single(getImages(StackName));
    else
        roiStack = single(getImages());
    end
catch
    roiStack = single(getImages());
    
end
% Get a Transmitted Image
% [nameTransFile,filePath] = uigetfile('*.tif','Select Transmitted Image for Overlay');
% imageTrans = imread([filePath,nameTransFile]);
clear check
for k = 1:size(files,1)
    current=files(k).name;
    if length(current)>=26
        check(k)=strcmp(current(end-25:end),'Transmitted Cell Image.tif');
    end
end
%%
loc=find(check);
if size(loc,1)==1
    imageTrans = imread(files(loc(1)).name);
else
    [nameTransFile,filePath] = uigetfile('*.tif','Select Transmitted Image for Overlay');
    imageTrans = imread([filePath,nameTransFile]);
end


% Get a Binary Image
% [nameBinary,filePath] = uigetfile('*.tif','Select Binary Image of Cell');
% imageBinary = imread([filePath,nameBinary]);
% [nameBinary2,filePath2] = uigetfile('*.tif','Select Binary Image of Shear Deformation');
% imageBinary2 = imread([filePath2,nameBinary2]);
% imageBinary2 = imresize(imageBinary2,[size(imageBinary,1) size(imageBinary,2)]);
clear check
for k = 1:length(files)
    current=files(k).name;
    if length(current)>=15
        check(k)=strcmp(current(end-14:end),'Binary Mask.tif');
    end
end
loc=find(check);
if size(loc,1)==1
    imageBinary= imread(files(loc(1)).name);
else
    [nameAreaFile,filePath] = uigetfile('*.tif','Select a Thresholded Image of the Cell Area');
    imageBinary = imread([filePath,nameAreaFile]);
end


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
%% Check to see if a previous row slope exists
clear files check
files = dir('*.mat'); %Check Directory for default filenames
filePath = strcat(cd,'\');
%Open Black Image
if size(files,1)>=1
    for k = 1:size(files,1)
        current=files(k).name;
        check(k)=strcmp(current(end-7:end),'rowV.mat');
    end
    loc=find(check);
else
    loc= zeros(1,1);
end
%%
if size(loc,2)<1
    w = msgbox('Select a location with low displacements and double-click to continue');
    waitfor(w);
    [~,sumBounds] = imcrop(sumImages);
    close
    
    sumBounds(1,3:4) = sumBounds(1,1:2) + sumBounds(1,3:4);
    sumBounds(1,5:6) = (sumBounds(1,1:2) + sumBounds(1,3:4))/2;
    sumBounds = sumBounds * xyScale;
end

boundsX = [0 size(imageBinary,2)*xyScale];
boundsY = [0 size(imageBinary,1)*xyScale];
%% Kilfoil Stack Filter
disp('Processing Images')
try
    load('processedStack.mat')
catch
    res=bpass3dMB(roiStack, [2 2 1], [7 7 12],[0 0]);
end
disp(['done Processing Images at '  num2str(toc) ' seconds'])
x = size(res,1);
y = size(res,2);
z = size(res,3);
%% Kilfoil Object Detection 3D
disp('Detecting 3D Centroids')
masscut = mean(prctile(roiStack(:,:,size(roiStack,3)),95));
r=feature3dMB(res, d , [d d 10], [x y z],[1 1 1],5,masscut,.3); %
r(:,1:2) = r(:,1:2)*xyScale;
r(:,3) = r(:,3)*zScale;
disp(['done Detecting 3D Centroids at ' num2str(toc) ' seconds'])
%%
save('rOrginal','r')
%% Build Planes Dot by Dot

clear rNbor
%Set a search window size and establish neighbors
radXY =  2.5; %microns
radZ = .3;
for i = 1:size(r,1)
    topX = r(i,1)+ radXY;
    botX = r(i,1)- radXY;
    topY = r(i,2)+ radXY;
    botY = r(i,2)- radXY;
    topZ = r(i,3)+ radZ;
    botZ = r(i,3)- radZ;
    rNbor(i,1:size(find(r(:,1)<topX & r(:,1)>botX & r(:,2)<topY & r(:,2)>botY& r(:,3)<topZ & r(:,3)>botZ))) = find(r(:,1)<topX & r(:,1)>botX & r(:,2)<topY & r(:,2)>botY& r(:,3)<topZ & r(:,3)>botZ);
end
%%

disp('Building Planes')
%Grow from starting point until no more plane members are found
clear planesTemp
working = 1;
searched = 1:1:size(r,1);
%start at first row in r
planesTemp(:,1) = rNbor(1,1:nnz(rNbor(1,:)));
planes = planesTemp;
j=1; %designates starting at plane 1
while working == 1  
    for i = 1:size(planes)
        if ismember(planes(i,1),searched) == 1
            clear new
            searched((planes(i,1)==searched)) = [];
            new(:,1) = rNbor(planes(i,1),1:nnz(rNbor(planes(i,1),:)));
            planesTemp = cat(1,planesTemp,new);
        end
    end
    sBefore = size(planes,1);
    planes = unique(cat(1,planes,planesTemp));
    sAfter = size(planes,1);
    if sBefore == sAfter
        planes2(1:size(planes,1),j) = planes(:,1);
        j=j+1;
        clear planes planesTemp
        for k = 1:size(r,1)
            if ismember(k,searched)==1
                planesTemp(:,1) = rNbor(k,1:nnz(rNbor(k,:)));
                planes = planesTemp;
                searched(searched==k) = [];
                break
            end
            if k == size(r,1)
                working = 0;
            end
        end
    end
end
disp(['done Building Planes at ' num2str(toc) ' seconds'])
%% View all detected planes
figure
hold on
for i = 1:size(planes2,2)
    scatter3(r(planes2(1:nnz(planes2(:,i)),i),1),r(planes2(1:nnz(planes2(:,i)),i),2),r(planes2(1:nnz(planes2(:,i)),i),3))
end
hold off
%% Filter planes with too few members (and update r)
clear planesFinal
j =1;
planesFinal = 0;
for i = 1:size(planes2,2)
    if nnz(planes2(:,i))>50
        planesFinal(1:nnz(planes2(:,i)),j) = planes2(1:nnz(planes2(:,i)),i);
        j=j+1;
    else
        for k = 1:nnz(planes2(:,i))
            r(planes2(k,i),:) =[];
            planes2((planes2>planes2(k,i))) = planes2((planes2>planes2(k,i)))-1;
            planesFinal((planesFinal>planes2(k,i))) = planesFinal((planesFinal>planes2(k,i)))-1;
        end
    end
end
%%
figure
hold on
for i = 1:size(planesFinal,2)
    scatter3(r(planesFinal(1:nnz(planesFinal(:,i)),i),1),r(planesFinal(1:nnz(planesFinal(:,i)),i),2),r(planesFinal(1:nnz(planesFinal(:,i)),i),3))
end
hold off


%% Dots in Cell Region
disp('Generating Row Slope')
%%
rNDC = zeros(1,1);
SE = strel('disk',50);
imageBinaryDilated = imerode(imageBinary,SE);
% imageBinaryHoriz = (sum(imageBinaryDilated,2))~=max(sum(imageBinaryDilated,2));
% imageBinaryVert  = sum(imageBinaryDilated,1)~=max(sum(imageBinaryDilated,1));
% imageBinaryBorder = ~kron(imageBinaryVert,imageBinaryHoriz);

% % imageBinaryComp = cat(3,imageBinary,imageBinaryDilated,imageBinaryBorder);
% imageBinaryDilated = uint8((uint8(imageBinary2).*imageBinaryDilated)+uint8(imerode((imageBinaryBorder==1),SE))>0);
% figure
% imshow(imageBinaryDilated,[])
for i = 1:size(r,1)
    %if it is under the cell
    if imageBinaryDilated(round(r(i,2)/xyScale),round(r(i,1)/xyScale))~=0
        rNDC = cat(1,rNDC,i);
    end
end
rNDC(1,:) = [];

%%
if size(loc,2)>0
    load('rowV.mat')
    disp('Found a previous row slope!')
else
    
    %% Dots in Cropped Region
    clear rND rNDB planesND planesNDNNZ planesNDIgnore planesDev neighbors
    %ND is non-deformed
    
    sumIndX  = (r(:,1)>sumBounds(1,1) & r(:,1)<sumBounds(1,3));
    sumIndY  = (r(:,2)>sumBounds(1,2) & r(:,2)<sumBounds(1,4));
    sumIndXY = sumIndX .* sumIndY;
    sumIndFinal = find(sumIndXY);
    rND = (r(sumIndFinal,1:end));
    ignoreCheck = 0;
    %Find plane with least fluctuations
    for i = 1:size(planesFinal,2)
        planesND(1:nnz(intersect(planesFinal(:,i),sumIndFinal)),i) = intersect(planesFinal(:,i),sumIndFinal);
    end
    for i = 1:size(planesND,2)
        planesNDNNZ(1,i) = nnz(planesND(:,i));
    end
    for i = 1:size(planesND,2)
        if nnz(planesND(:,i)) < .5*max(planesNDNNZ)
            ignoreCheck = 1;
            planesNDIgnore(1,i) = 1;
        end
    end
    if ignoreCheck == 1
        planesNDIgnore = (planesNDIgnore==0);
    else
        planesNDIgnore = ones(1,size(planesND,2));
    end
    
    %Determine best plane for an approximate non-deformed region
    for i = 1:size(planesFinal,2)
        planesDev(1,i) = std(r(planesND(1:nnz(planesND(:,i)),i),3));
    end
    planesDev = planesDev.*planesNDIgnore;
    planesBest = find(min(planesDev) == planesDev,1,'first');
    rNDB = r(planesND(1:nnz(planesND(:,planesBest)),planesBest),:);
    
    
    
    %Find Center dot (likely to have 4 equidistant neighbors)
    for i = 1:size(rNDB,1)
        differences(i,1) = rNDB(i,1)-sumBounds(1,5);
        differences(i,2) = rNDB(i,2)-sumBounds(1,6);
        differences(i,4) = sqrt(differences(i,1)^2 + differences(i,2)^2);
    end
    best = find(differences(:,4)==min(differences(:,4)));
    %%
    % figure
    % scatter3(rNDB(:,1),rNDB(:,2),rNDB(:,3))
    % hold on
    % scatter3(rNDB(best,1),rNDB(best,2),rNDB(best,3))
    % hold off
    %%
    k = 0;
    count = 0;
    while k == 0
        %Find the 4 'equidistant' neighbors
        clear differences sortedNew sortedOrig neighbors
        for i = 1:size(rNDB,1)
            differences(i,1:3) = rNDB(i,1:3)-rNDB(best,1:3);
            differences(i,4) = sqrt(differences(i,1)^2 + differences(i,2)^2);
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
    %%
    % figure
    % scatter3(rNDB(:,1),rNDB(:,2),rNDB(:,3))
    % hold on
    % scatter3(neighbors(:,1),neighbors(:,2),neighbors(:,3))
    % scatter3(rNDB(best,1),rNDB(best,2),rNDB(best,3))
    % hold off
    % % trckText = strcat('\leftarrow ',trckNum);
    % % text(lub(nghbrs(i,2),1)-(cntrPt(1,1)-fSizeXmin),lub(nghbrs(i,2),2)-(cntrPt(1,2)-fSizeYmin),lub(nghbrs(i,2),6),trckText,'Color','red')
    %
    %%
    m = 15;
    figure
    scatter3(rNDB(:,1),rNDB(:,2),rNDB(:,3))
    hold on
    scatter3(rNDB(sortedOrig(1:m,1),1),rNDB(sortedOrig(1:m,1),2),rNDB(sortedOrig(1:m,1),3))
    m=10;
    scatter3(rNDB(sortedOrig(1:m,1),1),rNDB(sortedOrig(1:m,1),2),rNDB(sortedOrig(1:m,1),3))
    m=5;
    scatter3(rNDB(sortedOrig(1:m,1),1),rNDB(sortedOrig(1:m,1),2),rNDB(sortedOrig(1:m,1),3))
    scatter(0,0,0)
     hold off
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
    %%
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
    %%
    if v1mean<v2mean
        rowV = (v1B-v1A)';
    else
        rowV = (v2B-v2A)';
    end
    
    rowV1 = (v1B-v1A)';
    rowV2 = (v2B-v2A)';
    %%
    figure
    hold on
    plot3([v1A(1) v1B(1)],[v1A(2) v1B(2)],[v1A(3) v1B(3)])
    plot3([v2A(1) v2B(1)],[v2A(2) v2B(2)],[v2A(3) v2B(3)])
    scatter3(0,0,0)
    hold off
    %% Display v1row and v2row
    
    figure
    hold on
    scatter3(rNDB(:,1),rNDB(:,2),rNDB(:,3))
    scatter3(rNDB(v1row(:,1),1),rNDB(v1row(:,1),2),rNDB(v1row(:,1),3))
    scatter3(rNDB(v2row(:,1),1),rNDB(v2row(:,1),2),rNDB(v2row(:,1),3))
    scatter3(0,0,0)
    hold off
    
end
disp(['done Generating Row Slope at ' num2str(toc) ' seconds'])

%%
disp('Building Rows')
[r,rows] = buildRows2(r,rowV,planesFinal);
disp(['done Building Rows at ' num2str(toc) ' seconds'])

%%
figure
hold on
for i = 1:size(rows,1)
    scatter3(r(rows(i,1:nnz(rows(i,:))),1),r(rows(i,1:nnz(rows(i,:))),2),r(rows(i,1:nnz(rows(i,:))),3))
end
%%
% Separate Rows by Plane
clear rowPlanes
for i = 1:size(planesFinal,2)
    for j = 1:size(rows,1)
        rowPlanes(j,1:size(intersect(rows(j,:),planesFinal(1:nnz(planesFinal(:,i)),i)),1),i) = intersect(rows(j,:),planesFinal(1:nnz(planesFinal(:,i)),i));
    end
end
%%
% Remake 'rows' variable with rows ordered by plane
clear newRows
nRS = find(rowPlanes(:,1,1)>0);
newRows(:,:) = rowPlanes(nRS,:,1);
rowPlanesIdx(1,1) = 1;
rowPlanesIdx(1,2) = size(nRS,1);
for i = 2:size(rowPlanes,3)
    clear currentPlanes
    nRS = find(rowPlanes(:,1,i)>0,1,'first');
    currentPlanes(:,:) =  rowPlanes((rowPlanes(:,1,i)>0),:,i);
    rowPlanesIdx(i,1) = size(newRows,1)+1;
    newRows = cat(1,newRows,currentPlanes);
    rowPlanesIdx(i,2) = size(newRows,1);
end
rows = newRows;

% Label Objects in 'r' with Their Row Number
for i = 1:size(rows,1)
    for j = 1:size(rows,2)
        if rows(i,j)>0
            r(rows(i,j),8) = i;
        end
    end
end

% Determine non-deformed members of each row
clear rowsNDCU rowsNDC
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
figure
hold on
for j = 1:size(rowPlanes,3)
    for i = 1:size(rowPlanes,1)
        n = nnz(rowPlanes(i,:,j));
        scatter3(r(rowPlanes(i,1:n,j),1),r(rowPlanes(i,1:n,j),2),r(rowPlanes(i,1:n,j),3))
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate plane distance from surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fS = open('fitSurface.mat');
try
    fitSurface = fS.fitSurface{3};
catch
    xtemp = [1;2;3;4;5]
    ytemp = [5;8;3;5;7]
    ztemp = size(res,3)*zScale*ones(5,1);
    fitSurface = fit([xtemp,ytemp],ztemp,'lowess','Span',0.1);
end
for j = 1:size(planesFinal,2)
    for i = 1:nnz(planesFinal(:,j))
        planesLoc(i,j) = (feval(fitSurface,r(planesFinal(i,j),1),r(planesFinal(i,j),2))) - r(planesFinal(i,j),3);
    end
end
planesLoc(planesLoc==0)=nan;
planesLoc2 = mean(planesLoc,'omitnan');

%%
% figure
% hold on
% for i = 1:size(planesFinal,2)
%     scatter3(r(planesFinal(1:nnz(planesFinal(:,i)),i),1),r(planesFinal(1:nnz(planesFinal(:,i)),i),2),r(planesFinal(1:nnz(planesFinal(:,i)),i),3))
% end
% plot(fitSurface)
% hold off
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2ND METHOD FOR MEASURING DISPLACEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Fitting Rows with Independent Slopes - Method 2')
rowFits = zeros(size(rows,1),3,2);
for i = 1:size(rows,1)
    if size(intersect(rows(i,:),rNDC),1)>1
        [A,B] = fitLine3D(r(intersect(rows(i,:),rNDC),1),r(intersect(rows(i,:),rNDC),2),r(intersect(rows(i,:),rNDC),3));
        rowFits(i,1:3,1) = A';
        rowFits(i,1:3,2) = B';
        
    elseif size(intersect(rows(i,:),rNDC),1)>0
        [xyzFinal,rowP] = transLine3D(rowV,r(intersect(rows(i,:),rNDC),1:3),imageSize);
        rowFits(i,1:3,1) = xyzFinal'-rowV;
        rowFits(i,1:3,2) = xyzFinal'+rowV;
    else
        [xyzFinal,rowP] = transLine3D(rowV,r(rows(i,1:nnz(rows(i,:))),1:3),imageSize);
        rowFits(i,1:3,1) = xyzFinal'-rowV;
        rowFits(i,1:3,2) = xyzFinal'+rowV;
    end
end
r2 = r;
%% Extend the fit Lines passed the edge of the dot area Method 2
clear rowFits2

ty = (boundsY(1,1) - rowFits(1,2,1))/(rowFits(1,2,2)-rowFits(1,2,1));
tx = (boundsX(1,1) - rowFits(1,1,1))/(rowFits(1,1,2)-rowFits(1,1,1));
for i = 1:size(rowFits,1)
    if abs(ty)<abs(tx)
        t = (boundsY(1,1) - rowFits(i,2,1))/(rowFits(i,2,2)-rowFits(i,2,1));
        rowFits2(i,1:3,1) = rowFits(i,1:3,1) + t*(rowFits(i,1:3,2)-rowFits(i,1:3,1));
        t = (boundsY(1,2) - rowFits(i,2,2))/(rowFits(i,2,1)-rowFits(i,2,2));
        rowFits2(i,1:3,2) = rowFits(i,1:3,2) + t*(rowFits(i,1:3,1)-rowFits(i,1:3,2));
    else
        t = (boundsX(1,1) - rowFits(i,1,1))/(rowFits(i,1,2)-rowFits(i,1,1));
        rowFits2(i,1:3,1) = rowFits(i,1:3,1) + t*(rowFits(i,1:3,2)-rowFits(i,1:3,1));
        t = (boundsX(1,2) - rowFits(i,1,2))/(rowFits(i,1,1)-rowFits(i,1,2));
        rowFits2(i,1:3,2) = rowFits(i,1:3,2) + t*(rowFits(i,1:3,1)-rowFits(i,1:3,2));
    end
end
rowFits2 = single(rowFits2);
%% Calculate Displacements from Fit Lines Method 2
rRef2a = zeros(size(r,1),3);
for i = 1:size(rowFits,1)
    if nnz(rowsNDCU(i,:))>0
        n = nnz(rows(i,:));
        temp = squeeze(rowFits2(i,:,:))';
        [xy,distance,t_a] = distance2curve(temp,r(rows(i,1:n),1:3));
        for j = 1:size(xy,1)
            rRef2a(rows(i,j),1:3) = xy(j,1:3);
        end
    else
        for j = 1:nnz(rows(i,:))
            rRef2a(rows(i,j),1:3) = r(rows(i,j),1:3);
        end
    end
end
%%
try
    load('XY Disp Data.mat')
catch
    book1 = -1;
end
 
book1 = single(book1);
book2 = single(book2);
clear book1N
r(:,9)=0;
frames = 1:1:size(book1,2);
for i = 1:size(book1,3)
book1(20,:,i) = frames';
book1(21,:,i) = i;
end
book1(book1==0) = NaN;
book1([1 2],:,:) = book1([1 2],:,:)*xyScale;
book1(20,:,:) = book1(20,:,:)*zScale;
book2(:,[1 2]) = book2(:,[1 2])*xyScale;

disp('Matching 3D detections to 2D-based pillars')
for i = 1:size(r,1)
    clear differences 
    differences = min(squeeze(sqrt((book1(1,:,:)-r(i,1)).^2+(book1(2,:,:)-r(i,2)).^2+(book1(20,:,:)-r(i,3)).^2)));
    
    if min(differences)<.5       
    r(i,9) = find(differences==min(differences));
    end

end
disp('done Matching 3D detections to 2D-based pillars')


%% Row Shift Correction
figure
hold on
for i=1:size(rows,1)
scatter3(r(rowsNDCU(i,(r(rowsNDCU(i,rowsNDCU(i,:)>0),9)>0),:),1),r(rowsNDCU(i,(r(rowsNDCU(i,rowsNDCU(i,:)>0),9)>0),:),2),r(rowsNDCU(i,(r(rowsNDCU(i,rowsNDCU(i,:)>0),9)>0),:),3))
end
%%

for i = 1:size(rows,1)
    rowsSCtest(i,1:length(r(rowsNDCU(i,(r(rowsNDCU(i,rowsNDCU(i,:)>0),9)>0),:),1))) = r(rowsNDCU(i,(r(rowsNDCU(i,rowsNDCU(i,:)>0),9)>0),:),1);
    rowsSC(i,1) = mean(r(rowsNDCU(i,(r(rowsNDCU(i,rowsNDCU(i,:)>0),9)>0),:),1)-book2(r(rowsNDCU(i,(r(rowsNDCU(i,rowsNDCU(i,:)>0),9)>0),:),9),1));
    rowsSC(i,3) = length(r(r(rowsNDCU(rowsNDCU(i,:)>0),9)>0,1)-book2(r(r(rowsNDCU(rowsNDCU(i,:)>0),9)>0,9),1));
    rowsSC(i,2) = mean(r(rowsNDCU(i,(r(rowsNDCU(i,rowsNDCU(i,:)>0),9)>0),:),2)-book2(r(rowsNDCU(i,(r(rowsNDCU(i,rowsNDCU(i,:)>0),9)>0),:),9),2));   
end

%store shift correction in r
for i = 1:size(r,1)
   r(i,10) = rowsSC(r(i,8),1);
    r(i,11) = rowsSC(r(i,8),2); 
end

%%
figure
hold on
for i=1:max(r(:,9))
    current = find(r(:,9)==i);
    if size(current,1)>0
        plot3(r(current,1),r(current,2),r(current,3))
    end
    
end
for i = 1:size(book1,3)
scatter3(book1(1,:,i),book1(2,:,i),book1(20,:,i),'.')
end
%% XY component of reference Method 2
for i = 1:size(r,1)
    
    temp = squeeze(rowFits2(r(i,8),1:2,:))';
    if r(i,9)>0
    [xy,distance,t_a] = distance2curve(temp,book2(r(i,9),1:2));
    rRef2b(i,1:2) = xy(1,1:2) + r(i,10:11);
    else
    rRef2b(i,1:2) = rRef2a(i,1:2);    
    end
end

%% Z component of reference Method 2    
for i = 1:size(r,1)
    differential = rowFits2(r(i,8),1:3,2)-rowFits2(r(i,8),1:3,1);
    differential2 = rowFits2(r(i,8),1:2,2) - rRef2b(i,1:2);
    differential3 = mean(differential(1,1:2)./differential2(1,1:2));
    rRef2b(i,3) = rowFits2(r(i,8),3,2)-differential(1,3)/differential3;
end
%%
figure
hold on
scatter3(rRef2b(:,1),rRef2b(:,2),rRef2b(:,3))
for i = 1:size(rowFits2,1)
    plot3([rowFits2(i,1,1) rowFits2(i,1,2)],[rowFits2(i,2,1) rowFits2(i,2,2)],[rowFits2(i,3,1) rowFits2(i,3,2)])
end

%%
rDisp2 = r(:,1:3)-rRef2b(:,1:3);
%% Displacement Statistics for Planes at least 4.5 microns from surface
% Find all dots farther than 4.5 microns from surface
clear rDisp2Mean
planesLocFilt = find(planesLoc2>4.5);
planesLocFiltList = planesFinal(:,planesLocFilt(1,1));
if size(planesLocFilt,2)>1
    for i = 2:size(planesLocFilt,2)
        planesLocFiltList = cat(1,planesLocFiltList,planesFinal(:,planesLocFilt(1,i)));
    end
end
planesLocFiltList(planesLocFiltList==0) = [];

%% Calculate Average Displacement per Row per Plane Method 2
[rDisp2Mean,rDisp2MeanPlanes,rDisp2StdPlanes,rDisp2MeanTotal,rDisp2StdTotal] = quickRowStats(rDisp2,planesFinal,rowPlanesIdx,rNDC,rowsNDCU,planesLocFiltList);
%% Filter out noise in Displacements Method 2
[noiseMean2,noiseStd2,noiseCutoff2,rDisp2Filt,rDisp2PF] = rowNoiseCalc(r,rDisp2,planesLocFiltList,rNDC);
disp('done Fitting Rows with Independent Slopes - Method 2')
%% Plot Questionable features
% figure
% %imshow(imageBinary,[])
% hold on
% scatter3(r(rDisp2PF(rDisp2PF(:,5)==1,1),1)/xyScale,r(rDisp2PF(rDisp2PF(:,5)==1,1),2)/xyScale,r(rDisp2PF(rDisp2PF(:,5)==1,1),3))
% scatter3(r(rDisp2PF(rDisp2PF(:,5)==0,1),1)/xyScale,r(rDisp2PF(rDisp2PF(:,5)==0,1),2)/xyScale,r(rDisp2PF(rDisp2PF(:,5)==0,1),3))
% scatter3(r(65,1),r(65,2),r(65,3))
% scatter3(0,0,0);

%% Scatter3/Plot3 of Dots/Fits Method 2
% figure
% hold on
% for j = 1:size(rowPlanes,3)
%     for i = 1:size(rowPlanes,1)
%         n = nnz(rowPlanes(i,:,j));
%         scatter3(r(rowPlanes(i,1:n,j),1),r(rowPlanes(i,1:n,j),2),r(rowPlanes(i,1:n,j),3))
%     end
% end
% 
% for i = 1:size(rowFits2,1)
%     plot3([rowFits2(i,1,1) rowFits2(i,1,2)],[rowFits2(i,2,1) rowFits2(i,2,2)],[rowFits2(i,3,1) rowFits2(i,3,2)])
% end
% scatter3(0,0,0)


%% Quiver Plot of Displacements Method 2
figure
%imshow(permute(res(:,:,1),[2,1,3]),[])
%hold on
quiver3(rRef2b(:,1),size(res,2)*xyScale-rRef2b(:,2),rRef2b(:,3),rDisp2(:,1),rDisp2(:,2)*-1,rDisp2(:,3),0)

%% Revisiting Method One reinforced by Average slopes
%Strategy: find closest and furthest points in row from one point in rowFit
%and use those two points to determine if a different fit using method 1 is
%necessary by querying whether either end point falls within the list of
%non-deformed dots under the cell.

% Determine what the ends of each row are
clear rowsEdgeDist rowsEnds
for i = 1:size(rows,1)
    clear distances currentRow currentMembers
    currentRow = rows(i,1:nnz(rows(i,:)));
    currentMembers = r(currentRow,1:3);
    for j = 1:size(currentMembers,1)
        rowsEdgeDist(i,j) = pdist([currentMembers(j,1:2);rowFits2(i,1:2,1)]);
    end
    rowsEnds(i,1) = rows(i,find(rowsEdgeDist(i,:) == min(rowsEdgeDist(i,1:nnz(rowsEdgeDist(i,:))))));
    rowsEnds(i,2) = rows(i,find(rowsEdgeDist(i,:) == max(rowsEdgeDist(i,1:nnz(rowsEdgeDist(i,:))))));
    rowsEnds(i,3) = pdist([r(rowsEnds(i,1),1:3);r(rowsEnds(i,2),1:3)]);
end

% Determine which rows are useful in calculating an average row slope
% Both ends must be outside of imageBinaryDilated's Cell and the line width
% must be at least 70 percent of maximum line width
for i = 1:size(rowPlanesIdx,1)
    rowPlanesMaxWidth(i,1) = max(rowsEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),3));
    rowsEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),5) = rowPlanesMaxWidth(i,1);
end
for i = 1:size(rowsEnds)
    if imageBinaryDilated(round(r(rowsEnds(i,1),2)/xyScale),round(r(rowsEnds(i,1),1)/xyScale))>0 && imageBinaryDilated(round(r(rowsEnds(i,2),2)/xyScale),round(r(rowsEnds(i,2),1)/xyScale))>0 && rowsEnds(i,3)>rowsEnds(i,5)*.7 %|| rowsEnds(i,?)>20
        rowsEnds(i,4) = 1;
    else
        rowsEnds(i,4) = 0;
    end
end

% Calculate an average row slope per plane
for i = 1:size(rowPlanesIdx,1)
    rowPlanesFits(i,1:3) = mean(rowFits2(find(rowsEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),4)==1),1:3,2)-rowFits2(find(rowsEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),4)==1),1:3,1));
    if size(find(rowsEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),4)==1),1) <1
        rowPlanesFits(i,1:3) = mean(rowFits2(find(rowsEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),4)==0),1:3,2)-rowFits2(find(rowsEnds(rowPlanesIdx(i,1):rowPlanesIdx(i,2),4)==0),1:3,1));
    end
end


% View ends and average fits
figure
hold on
scatter3(r(rowsEnds(:,1),1),r(rowsEnds(:,1),2),r(rowsEnds(:,1),3))
scatter3(r(rowsEnds(:,2),1),r(rowsEnds(:,2),2),r(rowsEnds(:,2),3))
scatter3(0,0,0)
plot3([0 rowPlanesFits(1,1)],[0 rowPlanesFits(1,2)],[0 rowPlanesFits(1,3)])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1ST METHOD FOR MEASURING DISPLACEMENTS V2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Fitting Rows with Row Slope - Method 1')
clear rowsLines
for i = 1:size(rowPlanesIdx,1)
    rowPlanesIdx2(i,1:pdist([rowPlanesIdx(i,1);rowPlanesIdx(i,2)])+1) = rowPlanesIdx(i,1):1:rowPlanesIdx(i,2);
end
%
for i = 1:size(rows,1)
    [rPIdxX rPIdxY] = find(rowPlanesIdx2 == i);
    if size(intersect(rows(i,:),rNDC),1)>0
        %[xyzFinal,rowP] = transLine3D(rowV,r(rowsNDCU(i,1:n),1:3),imageSize);
        [xyzFinal,rowP] = transLine3D(rowPlanesFits(rPIdxX,1:3),r(intersect(rows(i,:),rNDC),1:3),imageSize);
        rowsLines(1,1:3,i) = xyzFinal-rowPlanesFits(rPIdxX,1:3)';
        rowsLines(2,1:3,i) = xyzFinal+rowPlanesFits(rPIdxX,1:3)';
    else
        [xyzFinal,rowP] = transLine3D(rowPlanesFits(rPIdxX,1:3),r(rowsEnds(i,1:2),1:3),imageSize);
        rowsLines(1,1:3,i) = xyzFinal-rowPlanesFits(rPIdxX,1:3)';
        rowsLines(2,1:3,i) = xyzFinal+rowPlanesFits(rPIdxX,1:3)';
        i
    end
    
end

%% Extend the fit Lines passed the edge of the dot area - Method 1
% Can't believe this was a problem...

for i = 1:size(rowsLines,3)
    if mean(abs(rowPlanesFits(:,1)))<mean(abs(rowPlanesFits(:,2)))
        t = (boundsY(1,1) - rowsLines(1,2,i))/(rowsLines(2,2,i)-rowsLines(1,2,i));
        rowsLines2(1,1:3,i) = rowsLines(1,1:3,i) + t*(rowsLines(2,1:3,i)-rowsLines(1,1:3,i));
        t = (boundsY(1,2) - rowsLines(2,2,i))/(rowsLines(1,2,i)-rowsLines(2,2,i));
        rowsLines2(2,1:3,i) = rowsLines(2,1:3,i) + t*(rowsLines(1,1:3,i)-rowsLines(2,1:3,i));
    else
        t = (boundsX(1,1) - rowsLines(1,1,i))/(rowsLines(2,1,i)-rowsLines(1,1,i));
        rowsLines2(1,1:3,i) = rowsLines(1,1:3,i) + t*(rowsLines(2,1:3,i)-rowsLines(1,1:3,i));
        t = (boundsX(1,2) - rowsLines(2,1,i))/(rowsLines(1,1,i)-rowsLines(2,1,i));
        rowsLines2(2,1:3,i) = rowsLines(2,1:3,i) + t*(rowsLines(1,1:3,i)-rowsLines(2,1:3,i));
    end
end
%% Calculate Displacements from Fit Lines - Method 1
rRef1 = zeros(size(r,1),3);
for i = 1:size(rowsLines,3)
    if nnz(rowsNDCU(i,:))>0
        n = nnz(rows(i,:));
        [xy,distance,t_a] = distance2curve(rowsLines2(:,:,i),r(rows(i,1:n),1:3));
        for j = 1:size(xy,1)
            rRef1a(rows(i,j),1:3) = xy(j,1:3);
        end
    else
        for j = 1:nnz(rows(i,:))
            rRef1a(rows(i,j),1:3) = r(rows(i,j),1:3);
        end
    end
end

%%
clear differential differential2 differential3
for i = 1:size(r,1)
    
    temp = squeeze(rowsLines2(:,1:2,r(i,8)));
    if r(i,9)>0
    [xy,distance,t_a] = distance2curve(temp,book2(r(i,9),1:2));
    rRef1b(i,1:2) = xy(1,1:2) + r(i,10:11);
    else
    rRef1b(i,1:2) = rRef1a(i,1:2);    
    end
    
    differential = rowsLines2(2,:,r(i,8))-rowsLines2(1,:,r(i,8));
    differential2 = rowsLines2(2,1:2,r(i,8)) - rRef1b(i,1:2);
    differential3 = mean(differential(1,1:2)./differential2(1,1:2));
    rRef1b(i,3) = rowsLines2(2,3,r(i,8))-differential(1,3)/differential3;
    
end
%%
% figure
% hold on
% scatter3(rRef1b(:,1),(size(imageBinary,1)*xyScale)-rRef1b(:,2),rRef1b(:,3))
% for i = 1:size(rowsLines2,3)
%     plot3([rowsLines2(1,1,i) rowsLines2(2,1,i)],[(size(imageBinary,1)*xyScale)-rowsLines2(1,2,i) (size(imageBinary,1)*xyScale)-rowsLines2(2,2,i)],[rowsLines2(1,3,i) rowsLines2(2,3,i)])
% end

%%
rDisp1 = r(:,1:3)-rRef1b(:,1:3);
%% Calculate Average Displacement per Row per Plane Method 1
[rDisp1Mean,rDisp1MeanPlanes,rDisp1StdPlanes,rDisp1MeanTotal,rDisp1StdTotal] = quickRowStats(rDisp1,planesFinal,rowPlanesIdx,rNDC,rowsNDCU,planesLocFiltList);
%% Filter out noise in Displacements Method 1
[noiseMean1,noiseStd1,noiseCutoff1,rDisp1Filt,rDisp1PF] = rowNoiseCalc(r,rDisp1,planesLocFiltList,rNDC);

disp('done Fitting Rows with Row Slope - Method 1')
%% Scatter3/Plot3 of Dots/Fits Method 1
% figure
% hold on
% for j = 1:size(rowPlanes,3)
%     for i = 1:size(rowPlanes,1)
%         n = nnz(rowPlanes(i,:,j));
%         scatter3(r(rowPlanes(i,1:n,j),1),(size(imageBinary,1)*xyScale)-r(rowPlanes(i,1:n,j),2),r(rowPlanes(i,1:n,j),3))
%     end
% end
% for i = 1:size(rowsLines2,3)
%     plot3([rowsLines2(1,1,i) rowsLines2(2,1,i)],[(size(imageBinary,1)*xyScale)-rowsLines2(1,2,i) (size(imageBinary,1)*xyScale)-rowsLines2(2,2,i)],[rowsLines2(1,3,i) rowsLines2(2,3,i)])
% end
% scatter3(0,0,0)

%% Quiver Plot
% figure
% quiver3(rRef1b(:,1),rRef1b(:,2),rRef1b(:,3),rDisp1(:,1),rDisp1(:,2),rDisp1(:,3),0)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3RD METHOD FOR FITTING ROWS (METHOD 1 AND 2 COMBINED)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(rows,1)
    if rowsEnds(i,4) == 1
        rowFits3(i,1:3,1:2) = rowFits2(i,1:3,1:2);
    else
        rowFits3(i,1:3,1:2) = squeeze(rowsLines2(1:2,1:3,i))';
    end
end

%% Calculate Displacements from Fit Lines Method 3
rRef3a = zeros(size(r,1),3);
for i = 1:size(rowFits3,1)
    if nnz(rowsNDCU(i,:))>0
        n = nnz(rows(i,:));
        temp = squeeze(rowFits3(i,:,:))';
        [xy,distance,t_a] = distance2curve(temp,r(rows(i,1:n),1:3));
        for j = 1:size(xy,1)
            rRef3a(rows(i,j),1:3) = xy(j,1:3);
        end
    else
        for j = 1:nnz(rows(i,:))
            rRef3a(rows(i,j),1:3) = r(rows(i,j),1:3);
        end
    end
end
%%
for i = 1:size(r,1)
    
    temp = squeeze(rowFits3(r(i,8),1:2,:))';
    if r(i,9)>0
    [xy,distance,t_a] = distance2curve(temp,book2(r(i,9),1:2));
    rRef3b(i,1:2) = xy(1,1:2) + r(i,10:11);
    else
    rRef3b(i,1:2) = rRef3a(i,1:2);    
    end
    
    differential = rowFits3(r(i,8),1:3,2)-rowFits3(r(i,8),1:3,1);
    differential2 = rowFits3(r(i,8),1:2,2) - rRef3b(i,1:2);
    differential3 = mean(differential(1,1:2)./differential2(1,1:2));
    rRef3b(i,3) = rowFits3(r(i,8),3,2)-differential(1,3)/differential3;
    
end
%%
% figure
% hold on
% scatter3(rRef2b(:,1),rRef2b(:,2),rRef2b(:,3))
% for i = 1:size(rowFits3,1)
%     plot3([rowFits3(i,1,1) rowFits3(i,1,2)],[rowFits3(i,2,1) rowFits3(i,2,2)],[rowFits3(i,3,1) rowFits3(i,3,2)])
% end

%%
rDisp3 = r(:,1:3)-rRef3b(:,1:3);
%% Calculate Average Displacement per Row per Plane Method 3
[rDisp3Mean,rDisp3MeanPlanes,rDisp3StdPlanes,rDisp3MeanTotal,rDisp3StdTotal] = quickRowStats(rDisp3,planesFinal,rowPlanesIdx,rNDC,rowsNDCU,planesLocFiltList);
%% Filter out noise in Displacements Method 3
[noiseMean3,noiseStd3,noiseCutoff3,rDisp3Filt,rDisp3PF] = rowNoiseCalc(r,rDisp3,planesLocFiltList,rNDC);
disp('done Fitting Rows with Independent Slopes - Method 3')
%% Scatter3/Plot3 of Dots/Fits Method 3
% figure
% hold on
% for j = 1:size(rowPlanes,3)
%     for i = 1:size(rowPlanes,1)
%         n = nnz(rowPlanes(i,:,j));
%         
%         scatter3(r(rowPlanes(i,1:n,j),1),(size(imageBinary,1)*xyScale)-r(rowPlanes(i,1:n,j),2),r(rowPlanes(i,1:n,j),3))
%     end
%     
% end
% 
% for i = 1:size(rowFits3,1)
%     plot3([rowFits3(i,1,1) rowFits3(i,1,2)],[(size(imageBinary,1)*xyScale)-rowFits3(i,2,1) (size(imageBinary,1)*xyScale)-rowFits3(i,2,2)],[rowFits3(i,3,1) rowFits3(i,3,2)])
% end
% scatter3(0,0,0)
%%
% figure
% hold on
% scatter3(0,0,0)
% scatter3(x*xyScale,y*xyScale,z*zScale)
% for j = 1:size(rowPlanes,3)
%     if nnz(planes2(:,j))>5000
%         method3fit{j} = fit([r(planes2(1:nnz(planes2(:,j)),j),1),(size(imageBinary,1)*xyScale)-r(planes2(1:nnz(planes2(:,j)),j),2)],r(planes2(1:nnz(planes2(:,j)),j),3),'lowess','Span',.02);
%         plot(method3fit{j})
%     end
%     plot(fitSurface)
% end

%% Quiver Plot of Displacements Method 3
figure
hold on
for i = 1:size(planesFinal,2)
quiver3(rRef3b(planesFinal(1:nnz(planesFinal(:,i)),i),1),(size(imageBinary,1)*xyScale)-rRef3b(planesFinal(1:nnz(planesFinal(:,i)),i),2),rRef3b(planesFinal(1:nnz(planesFinal(:,i)),i),3),rDisp3(planesFinal(1:nnz(planesFinal(:,i)),i),1),rDisp3(planesFinal(1:nnz(planesFinal(:,i)),i),2)*-1,rDisp3(planesFinal(1:nnz(planesFinal(:,i)),i),3),0)
end
%% Create Histogram
close all
zdispdist = figure;
hold on
bins = [-400:20:400];
SigColor = [.1 .7 .7];
histogram(rDisp3(intersect(planesLocFiltList(:,1),rNDC),3)*1000,bins,'FaceColor',[.6 .6 .6],'Normalization','probability')
histmax = max(histcounts(rDisp3(intersect(planesLocFiltList(:,1),rNDC),3)*1000,bins,'Normalization','probability'))
p1 = plot([noiseCutoff3/2*1000 noiseCutoff3/2*1000],[0 round(histmax,2)+.01],'color',[.3 .3 .3],'linestyle','--','linewidth',1);
plot([(noiseCutoff3*-1)/2*1000 (noiseCutoff3*-1)/2*1000],[0 round(histmax,2)+.01],'color',[.3 .3 .3],'linestyle','--','linewidth',1)
p2 = plot([noiseCutoff3*1000 noiseCutoff3*1000],[0 round(histmax,2)+.01],'color',SigColor ,'linestyle','--','linewidth',1);
plot([noiseCutoff3*-1*1000 noiseCutoff3*-1*1000],[0 round(histmax,2)+.01],'color',SigColor ,'linestyle','--','linewidth',1)
set(gca,'fontsize',28)
xt = 'Z-Displacement (nm)';% input('enter the xaxis label','s');
yt = 'Probability'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = '\sigma'; %input('enter the legend','s');
le2 = '2*\sigma';
le3 = 'Cell Border';
xl = xlabel(xt);
yl = ylabel(yt);
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',28);
set(yl,'fontweight','bold','fontsize',28);
leg = legend([p1 p2],le,le2,'location','northwest');
leg.FontSize = 20;
axis([-400 400 0 round(histmax,2)+.01])

text((double(noiseCutoff3)*1.01*1000),(histmax*.5),strcat('2\sigma= ',num2str(round(noiseCutoff3*1000,0)),'nm'),'color',SigColor ,'fontsize',20)
text((double(noiseCutoff3)*1.01*1000),(histmax*.5-.01),'Noise Cutoff','color',SigColor ,'fontsize',18 )

title = '\Z-Displacement Histogram';
savefile = [filePath title];
export_fig(zdispdist,savefile,'-native');
%%
figure;
hold on
histogram(abs(rDisp3(intersect(planesLocFiltList(:,1),rNDC),3)),50)
histmax = max(histcounts(abs(rDisp3(:,3)),50));
plot([noiseCutoff3 noiseCutoff3],[0 histmax],'color',[.3 .3 .3],'linestyle','--','linewidth',1)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE COLOR MAPS OF DISPLACEMENTS IN Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorMapZ = single(brewermap(65536,'*PuBu'));
colorBar1 = single(zeros(500,25));
range = uint16(round(linspace(65536,1,500)'));
for i = 1:25
    colorBar1(1:500,i) = range;
end
colorBar2 = ind2rgb(colorBar1,colorMapZ);
for i = 1:10
    colorBar2((i*50)-3:(i*50),13:25,:) = 0;
end
colorBar2(1:3,13:25,:) = 0;

%Save Color Bar Image
close all
colorBarSave = figure;
hold on
imshow(colorBar2);
filePath=cd;
savefile = [filePath '\HeatMaps\ColorBarZv2.tif'];
export_fig(colorBarSave,savefile,'-native');

colorMapXY = single(brewermap(65536,'*Spectral'));
colorBar1 = single(zeros(500,25));
range = uint16(round(linspace(65536,1,500)'));
for i = 1:25
    colorBar1(1:500,i) = range;
end
colorBar2 = single(ind2rgb(colorBar1,colorMapXY));
for i = 1:10
    colorBar2((i*50)-3:(i*50),13:25,:) = 0;
end
colorBar2(1:3,13:25,:) = 0;

%Save Color Bar Image
close all
colorBarSave = figure;
hold on
imshow(colorBar2);
filePath=cd;
savefile = [filePath '\HeatMaps\ColorBarZv2.tif'];
export_fig(colorBarSave,savefile,'-native');


cutoff = 0.2; %microns
imageBinaryCombined = ((imageBinaryDilated==0));%+(imageBinary2==0))==0;
%%
%Determine which 'planes' in planesFinal should be the same plane
clear planesGroups
for i = 1:size(planesLoc2,2)
    clear differences
    differences = planesLoc2 - planesLoc2(1,i);
    planesGroups(i,1:size(find(abs(differences)<2),2)) = find(abs(differences)<2)';
end
planesGroups = unique(planesGroups,'rows');
[HeatMapN,vqN] = heatmapZ(r,rDisp3,planesFinal,planesGroups,imageBinaryCombined,xyScale,0,colorMapZ,colorMapXY,0);
[HeatMap,vq1] = heatmapZ(r,rDisp1Filt,planesFinal,planesGroups,imageBinaryCombined,xyScale,noiseCutoff1,colorMapZ,colorMapXY,1);
[HeatMap2,vq2] = heatmapZ(r,rDisp2Filt,planesFinal,planesGroups,imageBinaryCombined,xyScale,noiseCutoff2,colorMapZ,colorMapXY,2);
[HeatMap3,vq3] = heatmapZ(r,rDisp3Filt,planesFinal,planesGroups,imageBinaryCombined,xyScale,noiseCutoff3,colorMapZ,colorMapXY,3);


% Print Data to txt file
% Calculate Useful parameters
[vq1pos,vq1neg,rDisp1PosTotal, rDisp1PosMax, rDisp1NegTotal, rDisp1NegMax]  = vqStats(vq1,planesGroups);
[vq2pos,vq2neg,rDisp2PosTotal, rDisp2PosMax, rDisp2NegTotal, rDisp2NegMax]  = vqStats(vq2,planesGroups);
[vq3pos,vq3neg,rDisp3PosTotal, rDisp3PosMax, rDisp3NegTotal, rDisp3NegMax]  = vqStats(vq3,planesGroups);

% Write to file
cd HeatMaps
planesLocTxt = fopen('Average Z location of planes.txt','wt');
p1Format = 'Plane no. %1.0f is at %.2f microns from the surface \n';
p2Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 1 fit \n';
p3Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 2 fit \n';
p4Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 3 fit \n';
for i = 1:size(planesGroups,1)
    fprintf(planesLocTxt,p1Format,i,mean(planesLoc2(1,planesGroups(i,(planesGroups(i,:)>0)))));
    fprintf(planesLocTxt,p2Format,i,mean(rDisp1MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))));
    fprintf(planesLocTxt,p3Format,i,mean(rDisp2MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))));
    fprintf(planesLocTxt,p4Format,i,mean(rDisp3MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))));
end
fclose(planesLocTxt);
cd(filePath)

%Comma Delimiter Version

for j = 1:3
    cd HeatMaps
    planesLocComma = fopen(strcat('Planes Data ',num2str(j),'.txt'),'wt');
    p1Format = 'Plane\tDistance\tMean\tStd\tPositive Total\tPositive Max\tNegative Total\tNegative Max \n';
    p2Format = '%1.0f\t%.2f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f \n';
    fprintf(planesLocTxt,p1Format);
    for i = 1:size(planesGroups,1)
        fprintf(planesLocTxt,p2Format,i,mean(planesLoc2(1,planesGroups(i,(planesGroups(i,:)>0)))),mean(rDisp2MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))),mean(rDisp2StdPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))),rDisp2PosTotal(i,1),rDisp2PosMax(i,1),rDisp2NegTotal(i,1),rDisp2NegMax(i,1));
    end
    fclose(planesLocComma);
    cd(filePath)
end



%%
save 3DnormalData
%%
disp('Scipt has Completed')
% Save data for profile views
%%
filePath = cd;
folderName = 'Profile Data';
mkdir(filePath,folderName)
save('Profile Data\vqZ.mat','vqN','vq1','vq2','vq3','imageTrans','HeatMap','HeatMap2','HeatMap3')

%%
%save('rowV.mat','rowV')



%% Attempt to interpolate normal displacements in 3D
% rDispFilt = rDisp;
% rDispFilt((rDispFilt<.4 & rDispFilt>-.4))=0;
%
% meshRes = 2;
% %[xx,yy,zz] = meshgrid(1:meshRes:size(roiStack,1),1:meshRes:size(roiStack,2),1:1:size(roiStack,3));
% [xx,yy,zz] = meshgrid(1:meshRes:imageSize(1,1),1:meshRes:imageSize(2,1),min(r(:,3)):0.2:max(r(:,3)));
% xq = double([xx(:) yy(:) zz(:)]);
% vq = griddatan(rRef(:,1:3),rDispFilt(:,3),xq);
% vq = reshape(vq,size(xx));
% vq(isnan(vq)) = min(min(min(vq)));
%
% vq2 = vq;
% vq2 = vq2+abs(min(min(min(vq2))));
% vq2Scale = double(65000/max(max(max(vq2))));
% vq2 = uint16(vq2Scale*vq2);
% ShowStack(vq2,1,1)


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
toc