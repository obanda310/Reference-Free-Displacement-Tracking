%% Automated Script
clear all
close all
[ListName,ListPath] = uigetfile;
fullListName = strcat(ListPath,ListName);
load(strcat(ListPath,ListName),'dirList')
parts = strsplit(ListName, 'List');
prefix = parts{1};
set(0,'defaultfigurecolor',[1 1 1])
%% For updating paths (when files are moved)
%set a as new path root
%set b as old path root
ak = strfind(fullListName,ListName);

a = [fullListName(1:ak-1) ''];

Schar = char(dirList{:});
all_rows_same = all(diff(Schar, [], 1) == 0, 1);
common_cols = find(~all_rows_same, 1, 'first');
if isempty(common_cols)
    b = '?'
else
    b = dirList{1}(1:common_cols-1);
end

for i = 1:size(dirList,1)
    
    dirList{i,1} = strrep(dirList{i,1},b,a);
end
%%
save(fullListName,'dirList')

%% For clearing old files
for j = 1:2
    for i = 1:size(dirList,1)
        cleanOld(dirList{i,1})
        disp(num2str(i))
    end
end
%% TRACKING Script
clear failed
for i = 1:size(dirList,1)
    try
        tracking(dirList{i,1})
        disp(num2str(i))
    catch
        disp('Tracking Failed')
        i
    end
end
%% SHEAR deformation Scripts
clear failed
for i = 1:size(dirList,1)
    try
        dispShear(dirList{i,1})
        disp(num2str(i))
    catch
        disp('Shear Failed')
        i
    end
end

%% 3D deformation Scripts
clear failed3D
for i = 1:size(dirList,1)
    try
        disp3D(dirList{i,1})
        disp(num2str(i))
    catch
        disp('3D Failed')
        failed3D(i) = i;
        i
    end
end

%% Minimum Displacement Code
clear failedEmpty
for i = 1:size(dirList,1)
    try
        disp(num2str(i))
        dispEmpty(dirList{i,1})
        disp(num2str(i))
    catch
        disp('3D empty Failed')
        failedEmpty(i) = i;
        i
    end
end

%% Surface Displacement Code
clear failedSurface
for i = 1:size(dirList,1)
    try
        disp(num2str(i))
        dispSurface(dirList{i,1})
        disp(num2str(i))
    catch
        disp('3D surface Failed')
        failedSurface(i) = i;
        i
    end
end

%% For Running Stress Conversion
clear failedTractions TractionStats
failed = 0;
for i = 1:size(dirList,1)
    try
    disp(['Trying ' num2str(i)])
    [TractionStats(i,:)] = interpFinal3D_2(dirList{i,1});
    catch
        disp(['Conversion Failed at ' num2str(i)])
        failedTractions(i) = i;
    end

end
TractionStats(:,9:11) = TractionStats(:,4:6)*1000000000;
save([ListPath prefix 'TractionStats.mat'],'TractionStats','failed')
%% Visualizing Stats

results2 = TractionStats;
figure; scatter(results2(:,4),results2(:,5))
A = mean(results2(:,4)./results2(:,5),1,'omitnan')

    
%% Collect vqXY/vqZ totals for shear v normal graphs
%goodData = [1;2;3;5;6;14;17;23;31;34;43;44;47;49;50;48;11;12;19;25;28;36];
%goodData = [1:1:51]';

goodData = [1:1:size(dirList,1)]';
clear vqT
clear compZ
for m = 1:size(goodData,1)
    
    try
        clear planesLoc3
        disp(dirList{goodData(m,1),1})
        cd(dirList{goodData(m,1),1});
        load('3Ddata.mat')
        load('Profile Data\vqXY.mat');
        load('Profile Data\HeatMapXY.mat');
        %load('SurfaceData.mat','vqZt');
        load('Inputs2.mat','u')
        zTarget = 4;
        for j = 1:size(planesGroups,1)
            planesLoc3(j) = mean(planesLoc2(1,planesGroups(j,1:nnz(planesGroups(j,:)))));
        end
        zPlane = find(abs(planesLoc3-zTarget) == min(abs(planesLoc3-zTarget)),1,'first');
        vqXY(isnan(vqXY)) = 0;
        vqXY = imresize(vqXY,[size(vq3,1) size(vq3,2)]);
        m
        vqT(m,1) = sum(sum(vqXY));
        vq3(isnan(vq3)) = 0;
        vqT(m,2) = sum(sum(abs(u{1}{3}(:,:,end-1))))/(1*10^-6);
        vqT(m,3) = zPlane;
        vqT(m,4) = planesLoc2(zPlane);
        vqT(m,5) = sum(sum(image.Area ==0))*raw.dataKey(9,1)^2;
        vqT(m,6) = max(max(vqXY));
        vqT(m,7) = max(max(abs(vq3(:,:,zPlane))));
        
        compZ(1:500,1:500,m) = imresize(vq3(:,:,zPlane),[500 500],'nearest');
        compZ(501:1000,1:500,m) = imresize(vqZt(:,:),[500 500],'nearest');
        
    catch
        disp('failed')
    end
end

vqT(:,8) = vqT(:,1)./vqT(:,2);
vqT(:,9) = 1; vqT(:,1)>prctile(vqT(:,1),0);
vqTM = mean(vqT(vqT(:,9)==1,8));
vqTS = std(vqT(vqT(:,9)==1,8));
save([ListPath prefix 'ShearNormalStats.mat'],'vqT','vqTS','vqTM')
%% Compare Surface Z (dispSurface) deformation images to 3D-Z (disp3D) images
%need some way to qualitatively determine the accuracy of dispSurface
compZ(isnan(compZ)) = 0;
compZ = double(compZ);
compZ2 = (compZ+(-1*min(compZ(:))));
compZ2 = (compZ2/max(compZ2(:)))*256;
for i = 1:size(compZ,3)
    clear tempImage
    tempImage(:,:) = uint8(compZ2(:,:,i));
    
    imwrite(tempImage,'CompSurfaceTo3D.tif','WriteMode','append');
end
    


%% For Running Spherical Indentation Script
for i = 1:size(dirList,1)
    try
        SphericalIndent(dirList{i,1})
        disp(num2str(i))
    catch
        disp('Indentation Failed')
        i
    end
end


%% For Running Z Script
for i = 1:size(dirList,1)
    cd(dirList{i,1})
    %First check to see if the stack to analyze has been stored previously
    clear files
    files = dir('*.mat');
    if size(files,1)>=1
        for k = 1:length(files)
            current=files(k).name;
            if size(current,2)>12
                check(k)=strcmp(current(end-12:end),'StackName.mat');
            else
                check(k) = 0;
            end
        end
    end
    
    %If the Stack has not been assigned previously, manually assign it
    if size(find(check),2)==0
        [StackName,StackPath] = uigetfile('*.tif');
        save('StackName.mat','StackName')
    else
        load('StackName.mat','StackName')
        %if StackName == 0
            [StackName,StackPath] = uigetfile('*.tif');
            save('StackName.mat','StackName')
        %end
    end
    
end
%% Make sure every folder has a rowV file
[rowVName,rowVPath] = uigetfile;
load(strcat(rowVPath,rowVName),'rowV')
for i = 1:size(dirList,1)
    cd(dirList{i,1})
    save('rowV.mat','rowV')
end





%%
figure
scatter(vqT(:,1),vqT(:,2))
figure
scatter(vqT(:,5),vqT(:,1))
figure
scatter(vqT(:,5),vqT(:,2))
figure
scatter(vqT(:,6),vqT(:,7))
%% Imbalance Calculator
clear ftotals
for i = 1:size(dirList,1)
    clear totals
    totals = imbalance(dirList{i});
    ftotals(i,:) = totals;
    i
end
%%
ftotals(isnan(ftotals))=0;
(sum(ftotals(:,2)) + sum(ftotals(:,4)))/sum(ftotals(:,6))
%%

Ys(2,:) = ftotals(:,2)+ftotals(:,4);
Ys(3,:) = ftotals(:,1)+ftotals(:,3);
BarGrapher(Ys)



%% Calculate Overlap of Dilated Cell Area and Z-Deformation
%This code is to determine whether 8.125 micron dilation (50pixels@.1625)
%is sufficient to encompass all Z-deformation information
goodData = [1;2;3;5;6;14;17;23;31;34;43;44;47;49;50;48;11;12;19;25;28;36];
clear OverlapPct
for i = 1:size(goodData)
    OverlapPct(i) = nonDefOverlap(dirList{goodData(i)});
end


%% Max of each group

%% For calculating mean deviations in a dataset
clear StdVals
for i2 = 1:size(dirList,1)
    cd(dirList{i2,1})
    load('3Ddata.mat')
    xVals = shear.ltdX(:,shear.noCellTraj);
    xVals(xVals == 0) = NaN;
    
    yVals = shear.ltdY(:,shear.noCellTraj);
    yVals(yVals == 0) = NaN;
    
    xyVals = squeeze(shear.ltdXY(:,shear.noCellTraj));
    xyVals = cat(1,xyVals,-1*xyVals);
    xyVals(xyVals == 0) = NaN;
    
    StdVals(i2,1) = std(xVals(:),'omitnan')*1000;
    StdVals(i2,2) = std(yVals(:),'omitnan')*1000;
    StdVals(i2,3) = std(xyVals(:),'omitnan')*1000;
    StdVals(i2,4) = m3.noiseCutoff*1000;
    
    i2
    
end

xyzStats(1,1) = mean(StdVals(:,3)*2);
xyzStats(1,2) = std(StdVals(:,3)*2);

xyzStats(1,3) = mean(StdVals((StdVals(:,3)<(xyzStats(1,1)+2*xyzStats(1,2))),3)*2);
xyzStats(1,4) = std(StdVals((StdVals(:,3)<(xyzStats(1,1)+2*xyzStats(1,2))),3)*2);

xyzStats(2,1) = mean(StdVals(:,4));
xyzStats(2,2) = std(StdVals(:,4));

xyzStats(2,3) = mean(StdVals((StdVals(:,4)<(xyzStats(2,1)+2*xyzStats(2,2))),4));
xyzStats(2,4) = std(StdVals((StdVals(:,4)<(xyzStats(2,1)+2*xyzStats(2,2))),4));

