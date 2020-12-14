function dispEmpty(directory)
if nargin ==1
    cd(directory);
end
%%
clear all
close all
%%
tic
load('3Ddata.mat')

cols2 = r.col;
cols2(:,3) = 1:size(cols2,1);
cols2 = sortrows(cols2,[1 2]);
% cols2 Structure is [ShearPillar-ID, zframe, 3D-ID, dZ]
%%
for i = 1:size(cols2,1)-1
    if cols2(i,1) == cols2(i+1,1) && cols2(i,2)~=0
        cols2(i,4) = cols2(i+1,2)- cols2(i,2); %spacing between peaks
    else
        cols2(i,4) = NaN; %denotes the top 3D detection of a shear pillar
    end
end
coldist = mean(cols2(:,4),'omitnan'); %average spacing between peaks
if coldist == 0
    coldist = 10;
end
%%
clear emptyZ emptyZx emptyZy emptyZz emptyZid
optnew = optimset('TolX',1e-2);

%create filtered version of cols2 that only has pillars with multiple
%members
cols3 = cols2;
for i = 1:max(cols3(:,1))
    if size(find(cols3(:,1)==i),1)<=1
        cols3(cols3(:,1)==i,:) = [];
    end
end
cols3(cols3(:,1)==0,:) = [];
cols3u = unique(cols3(:,1));
progressbar('Finding Minimum Intensity Detections')
for i = 1:size(cols3,1)-1
    progressbar(i/size(cols3(:,1),1))
    if  abs((cols3(i+1,2)-cols3(i,2))-coldist)<coldist/2 && cols3(i+1,1)==cols3(i,1)         
        fun1 = fit([cols3(i,2):cols3(i+1,2)]',shear.Int(cols3(i,2):cols3(i+1,2),cols3(i,1)),'cubicinterp');
        emptyZz(i) = fminbnd(fun1,cols3(i,2),cols3(i+1,2),optnew);
        emptyZid(i) = i;
        if emptyZz(i)>0
            fun2 = fit([cols3(i,2):cols3(i+1,2)]',shear.rawX(cols3(i,2):cols3(i+1,2),cols3(i,1)),'cubicinterp');
            emptyZx(i) = feval(fun2,emptyZz(i));
            fun3 = fit([cols3(i,2):cols3(i+1,2)]',shear.rawY(cols3(i,2):cols3(i+1,2),cols3(i,1)),'cubicinterp');
            emptyZy(i) = feval(fun3,emptyZz(i));
        else
            emptyZx(i) = NaN;
            emptyZy(i) = NaN;
            emptyZz(i) = NaN;
        end
    end
end



try
emptyZ = cat(2,emptyZx',emptyZy',emptyZz',emptyZid');
emptyZ(emptyZ(:,1)==0,:)=[];
emptyZ(:,3) = emptyZ(:,3)*raw.dataKey(10,1);
catch
    emptyZ = zeros(1,4);
end
%%
%Try to build a layer at bottom of stack.
progressbar('Finding Bottom Minimum Intensity Detections')
count = 1;
for i = 1:max(cols2(:,1))
    progressbar(i/max(cols2(:,1)))
    if ismember(i,cols2(:,1))
    tempcol2 = cols2(find(cols2(:,1)==i,1,'first'),1:3);
    if tempcol2(1,2)>coldist*.75        
        fun1 = fit([1:tempcol2(1,2)]',shear.Int(1:tempcol2(1,2),i),'cubicinterp');
        emptyZb(count,3) = fminbnd(fun1,-1,(tempcol2(1,2)),optnew);
        emptyZb(count,4) = i;
        if emptyZb(count,3)>0
            fun2 = fit([1:tempcol2(1,2)]',shear.rawX(1:tempcol2(1,2),i),'cubicinterp');
            emptyZb(count,1) = feval(fun2,emptyZb(count,3));
            fun3 = fit([1:tempcol2(1,2)]',shear.rawY(1:tempcol2(1,2),i),'cubicinterp');
            emptyZb(count,2) = feval(fun3,emptyZb(count,3));
        else
            emptyZb(i,1:3) = NaN;
        end
        count = count+1;        
    end
    end
end

emptyZb(:,3) = emptyZb(:,3)*raw.dataKey(10,1);
%% Top Minimum Detections (does not work well, so not in use)
% 
% %Try to build a layer at top of stack.
% progressbar('Finding Top Minimum Intensity Detections')
% count = 1;
% for i = 1:max(cols2(:,1))
%     % try
%     progressbar(i/max(cols2(:,1)))
%     tempcol2 = cols2(find(cols2(:,1)==i,1,'last'),1:3);
%            
%         fun1 = fit([tempcol2(1,2):tempcol2(1,2)+ceil(coldist)]',shear.Int(tempcol2(1,2):tempcol2(1,2)+ceil(coldist),i),'cubicinterp');
%         emptyZc(count,3) = fminsearch(fun1,tempcol2(1,2)+ceil(coldist)/2);
%         emptyZc(count,4) = i;
%         if emptyZc(count,3)>0
%             fun2 = fit([tempcol2(1,2):tempcol2(1,2)+ceil(coldist)]',shear.rawX(tempcol2(1,2):tempcol2(1,2)+ceil(coldist),i),'cubicinterp');
%             emptyZc(count,1) = feval(fun2,emptyZc(count,3));
%             fun3 = fit([tempcol2(1,2):tempcol2(1,2)+ceil(coldist)]',shear.rawY(tempcol2(1,2):tempcol2(1,2)+ceil(coldist),i),'cubicinterp');
%             emptyZc(count,2) = feval(fun3,emptyZc(count,3));
%         else
%             emptyZc(i,1:3) = NaN;
%         end
%         count = count+1;        
%     
%     profileLines(1:1+ceil(coldist),i) = feval(fun1,tempcol2(1,2):tempcol2(1,2)+ceil(coldist));
%     % catch
%     % end
% end
% 
% figure
% hold on
% for i = 1:size(profileLines,2)
%     plot(1:1+ceil(coldist),profileLines(1:1+ceil(coldist),i))
% end
% xlim([0 1+ceil(coldist*0.7)])
% ylim([3 3.2])
% 
% emptyZc(:,3) = emptyZc(:,3)*raw.dataKey(10,1);

%% Combine Data
emptyZ = cat(1,emptyZb,emptyZ);
emptyZdel = find(abs(emptyZ(:,1))>max(size(image.Black)));
emptyZ(emptyZdel,:) = [];
emptyZdel = find(abs(emptyZ(:,1))==0);
emptyZ(emptyZdel,:) = [];
emptyZdel = find((emptyZ(:,1))<0);
emptyZ(emptyZdel,:) = [];
emptyZ(isnan(emptyZ(:,1)),:) = [];


%%
figure
scatter3(emptyZ(:,1), emptyZ(:,2), emptyZ(:,3),'b')
xlim([0 max(r.r(:,1))+5])
ylim([0 max(r.r(:,2))+5])
zlim([0 size(shear.rawX,1)*raw.dataKey(10,1)])
hold on
scatter3(r.r(:,1), r.r(:,2), r.r(:,3),'r')

%%
figure
plot([12:24]',shear.Int(12:24,60))
hold on
fun1 = fit([12:24]',shear.Int(12:24,60),'cubicinterp');
plot(fun1)
scatter(fminsearch(fun1,18),feval(fun1,17.6209))
%% Run disp3D protocol to collect zero positions and measure displacements
%Build Raw 3D Object
erraw = RawData3D(res,raw,emptyZ(:,1:3));
erraw = TranscribeR(erraw);

%Build Planes
planeer = PlanesData(erraw);
planeer = nborsPlanesF(planeer,erraw,radXY,radZ);
planeer = growPlanes(planeer,erraw);
[planeer,er] = cleanPlanes(planeer,erraw);
er = RawData3D(res,raw,er);
er = TranscribeR(er);


%Identify Non-Deformed Detections
er = regionCheck(er,image.ADil,raw);

%Build Rows
[er,errows] = buildRows2(er,rowV,planeer.final);

%% Final Cleaning Step - Removing all single-member row objects 
% This will help with processing time on row fits and should clean up some
% noise in poor quality datasets.
clear rowsNum rowspurge
for i = 1:size(errows,1)
rowsNum(i) = nnz(errows(i,:));
end
rowspurge = errows((rowsNum==1)',1);

rTemp = er.r;
rTemp(rowspurge,:) = [];

% Refresh raw data variable and planes data
er = RawData3D(res,raw,rTemp);
er = TranscribeR(er);

clear plane
planeer = PlanesData(er);
planeer = nborsPlanesF(planeer,er,radXY,radZ);
%Grow from starting point until no more plane members are found
planeer = growPlanes(planeer,er);
[planeer,er] = cleanPlanes(planeer,er);
er = RawData3D(res,raw,er);
er = TranscribeR(er);
er = regionCheck(er,image.ADil,raw);

disp('Rebuilding Rows')
[er,errows] = buildRows2(er,rowV,planeer.final);
disp(['done rebuilding Rows at ' num2str(toc) ' seconds'])



%Format Rows
[errows,errowPlanes,errowPlanesIdx,errowsNDCU,er]=formatRows(errows,planeer,er);

%Calculate plane distance from surface
[erplanesLoc2,erplanesLocFiltList] = placePlanes(er,planeer,Surface2);

%Match to 2D Pillars
[er] = match2D(er,raw,shear,errowsNDCU);

%1st Fit Method
em1 = DispData3D;
em1 = method1fit(em1,er,errows,rowV);
em1 = calcDisp(em1,er,errows,errowsNDCU); %based on fits of non-deformed markers
em1 = calcDispSC(em1,er,shear); %includes shift correction in reference approximation
em1.disp = er.r(:,1:3)-em1.refSC(:,1:3);
em1 = dispStats(em1,planeer,errowPlanesIdx,er,errowsNDCU,erplanesLocFiltList);
em1 = dispNoise(em1,er,erplanesLocFiltList,image,shear,raw.dataKey(9,1));

%2nd Fit Method
em2 = DispData3D;
em2 = method2fit(em2,em1,er,raw,image.ADil,errows,errowPlanesIdx);
em2 = calcDisp(em2,er,errows,errowsNDCU); %based on fits of non-deformed markers
em2 = calcDispSC(em2,er,shear); %includes shift correction in reference approximation
em2.disp = er.r(:,1:3)-em2.refSC(:,1:3);
em2 = dispStats(em2,planeer,errowPlanesIdx,er,errowsNDCU,erplanesLocFiltList);
em2 = dispNoise(em2,er,erplanesLocFiltList,image,shear,raw.dataKey(9,1));

%3rd Fit Method
em3 = DispData3D;
em3 = method3fit(em3,em1,em2,errows);
em3 = calcDisp(em3,er,errows,errowsNDCU); %based on fits of non-deformed markers
em3 = calcDispSC(em3,er,shear); %includes shift correction in reference approximation
em3.disp = er.r(:,1:3)-em3.refSC(:,1:3);
em3 = dispStats(em3,planeer,errowPlanesIdx,er,errowsNDCU,erplanesLocFiltList);
em3 = dispNoise(em3,er,erplanesLocFiltList,image,shear,raw.dataKey(9,1));
%%
% HeatMaps
clear erplanesGroups
for i = 1:size(erplanesLoc2,2)
    clear differences
    differences = erplanesLoc2 - erplanesLoc2(1,i);
    erplanesGroups(i,1:size(find(abs(differences)<2),2)) = find(abs(differences)<2)';
end
erplanesGroups = unique(erplanesGroups,'rows');

[HeatMapN,vqN] = heatmapZ(er.r,em3.disp,planeer.final,erplanesGroups,imageBinaryCombined,image.Borders,raw.dataKey(9,1),0,colorMapZ,colorMapXY,0);
%[HeatMap3,vq3] = heatmapZ(er.r,em3.disp,planeer.final,erplanesGroups,imageBinaryCombined,image.Borders,raw.dataKey(9,1),em3.noiseCutoff,colorMapZ,colorMapXY,3);
[HeatMap3,vq3] = heatmapZ(er.r,em3.disp,planeer.final,erplanesGroups,imageBinaryCombined,image.Borders,raw.dataKey(9,1),em3.noiseCutoff,colorMapZ,colorMapXY,3);

%%
save('empties.mat','emptyZ','er','em3','planeer')
%%
toc