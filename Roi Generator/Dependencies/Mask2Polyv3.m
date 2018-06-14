function [poly1, poly2, poly3] = Mask2Polyv3(raw,options)

xsize = size(raw,2);
ysize = size(raw,1);
LSS =1; % Change this value to change how Version 2 (old version) works.
%%
% Raw = binary mask to be converted to polygons
% LSS = Final edge size upon scale down of raw image. Higher number means
% higher resolution, but increased number of vertices per polygon.


% raw(:,1:EDGE) = 0;
% raw(1:EDGE,:) = 0;
% raw(:,end-EDGE:end) = 0;
% raw(end-EDGE:end,:) = 0;


%% Identify closed regions in image
rawN = raw==0;
%Remove small closed regions
rawN2 = bwareaopen(rawN, round(xsize/50)^2, 4);
%replace the original image 'raw'
raw = rawN2==0;


%% Remove isolated features up to 10pixels area.
if options(1,1) == 1
raw = bwareaopen(raw, 10, 4);
end

%% Open closed regions
% This is difficult to code well...could use some work.

%% Divide the image to break open large closed spaces

div = 10;
IMsize = ysize;
IMdiv = round(IMsize/div);
DivMask= ones(ysize,xsize);

Dz= zeros(ysize,xsize);
D = Dz;
% The image 'D' will contain only the edges created by dividing the raw
% image (used in v2)
if options(4,1) == 1
for i =1:div-1
    DivMask(IMdiv*i,:) = 0;
    %A(:,IMdiv*i) = 0;
    D(IMdiv*i-1,:) = 1;
    D(IMdiv*i+1,:) = 1;    
end
end

%% Identify features that do not need to be broken up.
smallraw = xor(bwareaopen(raw,round((xsize/13)^2),4),raw);
DivMask2 = smallraw | DivMask;

%% Update 'raw'
IM2 = (raw.*DivMask2); %Update 'raw' image by convolution with DivMask
poly1 = bwboundaries(IM2,8); %Find boundaries to start building polygons

%% %%%%%%%%%%%%%%%%%%%%%%%%%%VERSION 2 Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 2 used a method where the user had the option to scale down the
% image to generate a simplified mask and therefore simpler regions. It may
% still have use, but is currently not used.
if options(3,1) ==1
% %%Identify 'corners' in a simplified image
points = detectHarrisFeatures(imresize(raw,[ysize*LSS xsize*LSS]));
points2 = points;
points2.Location = round(points.Location/LSS);

% %% Find corners created by dividing the image
IM3 = IM2.*D; %Get synthetic edges in 'raw'
D2 = Dz;
%Mark in D2 anything that doesn't look like the end of a line in IM3
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

% %% Create an image of the outlines
D5 = Dz;
for i = 1:size(poly1,1)
    thisbound = poly1{i,1};
    for j = 1:size(thisbound,1)
        D5(thisbound(j,1),thisbound(j,2)) = 1;
    end
end

% %% Create an image of the outlines minus synthetic edges
D6 = D5.*(IM3==0);
% figure
% imshow(D6)

% %% For each 'corner' found in the simplied image (i.e. corners in
% points2), identify a corresponding pixel in the outlines in 'D6' and
% record that pixel in 'ptsFinal'
[idxY,idxX] = find(D6);
for i = 1:size(points2.Location,1)
    clear idxD
    idxD = (idxX-points2.Location(i,1)).^2+(idxY-points2.Location(i,2)).^2;
    idxM = find(idxD == min(idxD),1,'first');
    ptsFinal(i,2)=idxX(idxM);
    ptsFinal(i,1)=idxY(idxM);
    
end

% %% For each point in ptsFinal, identify the corresponding point in poly1
for j = 1:size(poly1,1)
    ptsInt = intersect(poly1{j,1}(:,1:2),ptsFinal,'rows');
    for i = 1:size(ptsInt,1)
        [~,idx] = ismember(ptsInt(i,1:2),poly1{j,1}(:,1:2),'rows');
        poly1{j,1}(idx(1,1),3) = 1;
    end
end

% %% Collect corners caused by dividing the image
% Identify locations in image correspond to corners
IM4 = IM3.*(D2==0);
[corners(:,1),corners(:,2)] = find(IM4);
% figure
% imshow(IM4,[])

% Identify which pixels in poly1 correspond to corners
for j = 1:size(poly1,1)
    ptsInt = intersect(poly1{j,1}(:,1:2),corners,'rows');
    for i = 1:size(ptsInt,1)
        [~,idx] = ismember(ptsInt(i,1:2),poly1{j,1}(:,1:2),'rows');
        poly1{j,1}(idx(1,1),3) = 1;
    end
end

% %% Keep only identified features in 'poly1' outlines.
for i = 1:size(poly1,1)
    try
        idx = find(poly1{i,1}(:,3));
        poly2{i,1} = poly1{i,1}(idx,1:2);
    end
end

% %% Begin converting data into format for creating .Regions file
for i = 1:size(poly2,1)
    spoly2(i,1) = size(poly2{i,1},1);
end
else
    poly2 = 'Change the option to run Mask2Poly version 2 if you want this output';
    corners = [0,0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%VERSION 2 End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% VERSION 3
%Builds off some ideas of v2. Attempts to get rid of straight lines instead
%of relying on corner finding algorithms, which don't seem to be exclusive
%enough

for i2 = 1:size(poly1,1)
    %% Start by setting pT to the boundary of a polygon "i2" from poly
    clear pT3 pT4 pT5 pT6
    pT3 = zeros(size(poly1{i2,1},1),6);
    pT3(:,1:2) = poly1{i2,1}(:,1:2);
    ptsInt = intersect(pT3(:,1:2),corners,'rows');
    for i = 1:size(ptsInt,1)
        [~,idx] = ismember(ptsInt(i,1:2),pT3(:,1:2),'rows');
        %[~,idx] = ismember(,);
        pT3(idx(1,1),6) = 1;
    end
    
    
    %% Identify Pixels in Vertical or Horizontal Lines
    % Also making sure not to mark corners
    for i = 1:size(pT3,1)-2
        
        
        if ((pT3(i,1) == pT3(i+1,1) && pT3(i+1,2) ~= pT3(i+2,2)) || (pT3(i,2) == pT3(i+1,2) && pT3(i+1,1) ~= pT3(i+2,1))) && pT3(i,6)==0
            pT3(i+1,3) = 1;
        end
        if pT3(i,1) == pT3(i+2,1) && pT3(i,2) == pT3(i+2,2)
            pT3(i+1,6) = 1;
            pT3(i+1,3) = 0;
        end
    end
    
    % Remove Pixels in Vertical or Horizontal Lines
    index = find(pT3(:,3)==0);
    pT4 = pT3(index,1:6);
    pT4(:,3:5) = 0;
    
    %% Identify slope from one pixel to next
    for i = 1:size(pT4,1)
        if i == size(pT4,1)
            pT4(i,3) = pT4(i,1)-pT4(1,1);
            pT4(i,4) = pT4(i,2)-pT4(1,2);
        else
            pT4(i,3) = pT4(i,1)-pT4(i+1,1);
            pT4(i,4) = pT4(i,2)-pT4(i+1,2);
        end
    end
    % Identify consecutive pixels with same slope
    for i = 1:size(pT4,1)-1
        if pT4(i,3) == pT4(i+1,3) && pT4(i,4) == pT4(i+1,4) && pT4(i+1,6) == 0
            pT4(i+1,5) = 1;
        end
    end
    % Remove consecutive pixels with same slope
    
    index = find(pT4(:,5)==0);
    pT5 = pT4(index,1:6);
    pT5(:,5)=0;
    %% Identify consecutive pairs of pixels with same slope
    for i = 1:size(pT5,1)-3
        if sum(pT5(i:i+1,3) == pT5(i+2:i+3,3))==2 && sum(pT5(i:i+1,4) == pT5(i+2:i+3,4))==2 && pT5(i+2,6) ==0 && pT5(i+3,6) ==0
            pT5(i+2:i+3,5) = 1;
        end
    end
    % Remove consecutive pairs of pixels with same slope
    
    index = find(pT5(:,5)==0);
    pT6 = pT5(index,1:6);
    pT6(:,1:2) = round(pT6(:,1:2));
    poly3{i2,1}(:,1) = pT6(:,2);
    poly3{i2,1}(:,2) = pT6(:,1);
end

%% View all polylines
% figure
% imshow(raw)
% hold on
% for i = 1:size(poly3,1)
%     try
% plot([poly3{i,1}(:,1);poly3{i,1}(1,1)],[poly3{i,1}(:,2);poly3{i,1}(1,2)])
% scatter([poly3{i,1}(:,1);poly3{i,1}(1,1)],[poly3{i,1}(:,2);poly3{i,1}(1,2)])
%     end
% end

%% Calculate total number of vertices
%if no regions, skip this
try
for i = 1:size(poly3,1)
    elemsizes(i,1) = size(poly3{i,1},1);
end
numVertices = sum(elemsizes,1);
catch
    numVertices = 0;
    poly1 = 0;
    poly2 = 0;
    poly3 = 0;
    disp('Warning! No vertices on current image! If that was intended, ignore this!')
end
