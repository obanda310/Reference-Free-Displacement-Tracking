function [polygons,C] = Mask2Polyv2(raw,LSS)
%%
% Raw = binary mask to be converted to polygons
% LSS = Final edge size upon scale down of raw image. Higher number means
% higher resolution, but increased number of vertices per polygon.


xsize = size(raw,2);
ysize = size(raw,1);
%Make sure there are no features within 'EDGE' pixels of border
EDGE = round(4000/LSS);
raw(:,1:EDGE) = 0;
raw(1:EDGE,:) = 0;
raw(:,end-EDGE:end) = 0;
raw(end-EDGE:end,:) = 0;

%%
points = detectHarrisFeatures(imresize(raw,[LSS LSS]));
points2 = points;
points2.Location = round(points.Location* (xsize/LSS));

% 
% figure
% imshow(raw)
% hold on
% plot(points2)

div = 10;
IMsize = xsize;
IMdiv = round(IMsize/div);
A2= ones(IMsize,IMsize);
Dz= zeros(IMsize,IMsize);
D = Dz;
D2 = Dz;

E= zeros(IMsize+2,IMsize+2);
for i =1:div-1
A2(IMdiv*i,:) = 0;   
%A(:,IMdiv*i) = 0;
D(IMdiv*i-1,:) = 1;
D(IMdiv*i+1,:) = 1;

end
IM2 = (raw.*A2);
IM2L = bwlabel(IM2);
IM3 = IM2.*D;
for i =1:div-1
for j = 2:xsize-1
    if IM3(IMdiv*i-1,j-1)==1 && IM3(IMdiv*i-1,j+1)==1 
        D2(IMdiv*i-1,j)=1;
    end
    if IM3(IMdiv*i+1,j-1)==1 && IM3(IMdiv*i+1,j+1)==1 
        D2(IMdiv*i+1,j)=1;
    end
end
end
IM4 = IM3.*(D2==0);

E(2:xsize+1,2:xsize+1) = IM3;

points3 = detectHarrisFeatures(IM3);
C = bwboundaries(IM2,8);
% 
% figure
% imshow(IM2)
% hold on
% plot(points2)
% 
% figure
% imshow(IM4,[])
% hold on 
% plot(points3)

%%
D5 = Dz;
for i = 1:size(C,1)
   thisbound = C{i,1};
   for j = 1:size(thisbound,1)
   D5(thisbound(j,1),thisbound(j,2)) = 1;
   end
end

D6 = D5.*(IM3==0);
% figure
% imshow(D6)


Dz= zeros(IMsize,IMsize);
se = strel('disk',10);
clear ptsFinal
[idxY,idxX] = find(D6);
for i = 1:size(points2.Location,1)
clear idxD
idxD = (idxX-points2.Location(i,1)).^2+(idxY-points2.Location(i,2)).^2;
idxM = find(idxD == min(idxD),1,'first');
ptsFinal(i,2)=idxX(idxM);
ptsFinal(i,1)=idxY(idxM);

end
%%
% figure
% imshow(raw)
% hold on
% scatter(ptsFinal(:,2),ptsFinal(:,1))
% scatter(C{1,1}(:,2),C{1,1}(:,1))

%%
for j = 1:size(C,1)
    ptsInt = intersect(C{j,1}(:,1:2),ptsFinal,'rows');
    for i = 1:size(ptsInt,1)
       % idx = find((C{j,1}(:,1:2))==(ptsInt(i,1:2)),1);
       [~,idx] = ismember(ptsInt(i,1:2),C{j,1}(:,1:2),'rows'); 
        C{j,1}(idx(1,1),3) = 1;
    end
end

%%
clear corners
[corners(:,1),corners(:,2)] = find(IM4);

for j = 1:size(C,1)
    ptsInt = intersect(C{j,1}(:,1:2),corners,'rows');
    for i = 1:size(ptsInt,1)
        [~,idx] = ismember(ptsInt(i,1:2),C{j,1}(:,1:2),'rows');
        %idx = find((C{j,1}(:,1:2))==(ptsInt(i,1:2)),1);
       %[~,idx] = ismember(,); 
        C{j,1}(idx(1,1),3) = 1;
    end
end
%%
for i = 1:size(C,1)
    try
    idx = find(C{i,1}(:,3));
    C2{i,1} = C{i,1}(idx,1:2);
    end
end


%%
figure
imshow(raw)
hold on
for i = 1:size(C2,1)
    try
plot(C2{i,1}(:,2),C2{i,1}(:,1))
    end
end

%% Begin converting data into format for creating .Regions file
for i = 1:size(C2,1)
    sC2(i,1) = size(C2{i,1},1);
end
NVert = sum(sC2);
numElements = size(C2,1);
polygons = C2;
end
