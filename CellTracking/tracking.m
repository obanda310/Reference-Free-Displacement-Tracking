%% This needs its own section
clear;
close all;
clc;

%% Get Images and Metadata
experiment = Experiment;
images = experiment.images;
roiCell = experiment.cellImg;
meta = experiment.metadata;
noImgs = size(images,3);
% Scaling is the same in X and Y; convert from meters to microns
pixelSize = meta.scalingX*1000000;

%% Final Pre-processing Before Finding Local Maxima
clear roiCell;
[roiImgs,roiMasks,roiCell,roiBounds,bkImg] = experiment.cropImgs;
scaleFactor = pixelSize/0.0825
%% Finding 3D Local Maxima
% Create a gaussian filtered version of original to decrease false local
% maxima
d = 3;
sig1 = 1/(1+sqrt(2))*d;
sig2 = sqrt(2) * sig1;

% Multiply the gaussian image by the mask image to isolate regions of
% interest
ppImages7 = double(imgaussfilt(roiImgs,sig2));
ppImages8 = roiMasks.*ppImages7;
ShowStack(ppImages8,experiment.centroids2d)

% Find local maxima in 3D (pixel resolution)
localMaxima3D = imregionalmax(ppImages8);

%% 2D maxima approach
% Taking a break from working with 3D maxima because too many data points
% are lost in the process, and it is seeming like it will not be a good way
% to eventually identify ellipsoids and their strain/displacement. Will now
% attempt to create ellipsoids by identifying all local 2D maxima belonging
% to a single pillar, and tracing the the major axis of individual
% ellipsoids through local maxima.
clear rM cM sM rsM csM ppImages10 subpixMaxima
clear tempInd tempInd2 tempInd3 tempInd4
close all
% Find local maxima in 2D (pixel resolution)
localMaxima2D = zeros(size(ppImages8));
for i = 1:noImgs
    thisMax2D = imregionalmax(ppImages8(:,:,i));    
    % In the event that a frame is empty, the local maxima are the entire
    % image (0's), this if statement removes these maxima.
    if thisMax2D == ones(size(thisMax2D,1),size(thisMax2D,2))
        thisMax2D = zeros(size(thisMax2D,1),size(thisMax2D,2));
    end
    localMaxima2D(:,:,i) = thisMax2D;
end
[rM,cM,sM] = ind2sub(size(localMaxima2D),find(localMaxima2D == 1));

% Separate maxima by z (frame) and determine indices in rM, cM, sM
for i = 1:noImgs
    if min(find(sM == i)) > 0
        twoDimMaxInd(i,1) = min(find(sM == i));
        twoDimMaxInd(i,2) = max(find(sM == i));
    end
end

% Use indices to create book of subpixel maxima for use later in linking
% maxima to pillars
for i = 1:noImgs
    if min(find(sM == i)) > 0
        [rsM,csM] = subpix2d(rM(twoDimMaxInd(i,1):twoDimMaxInd(i,2)),cM(twoDimMaxInd(i,1):twoDimMaxInd(i,2)),double(ppImages7(:,:,i)));
        subpixMaxima(1:size(rsM,2),1,i) = rsM(1,:);
        subpixMaxima(1:size(rsM,2),2,i) = csM(1,:);
        subpixMaxima(1:size(rsM,2),3,i) = i;
    end
end

%clear out-of-bounds results from subpix2d

%-greater than x image size
[tempInd, tempInd2] = find(subpixMaxima(:,1,:) > size(ppImages8,1));
subpixMaxima(tempInd,1:3,tempInd2) = 0;
%-greater than y image size
[tempInd, tempInd2] = find(subpixMaxima(:,2,:) > size(ppImages8,2));
subpixMaxima(tempInd,1:3,tempInd2) = 0;
%-smaller than 0 in x
[tempInd, tempInd2] = find(subpixMaxima(:,1,:) < 0);
subpixMaxima(tempInd,1:3,tempInd2) = 0;
%-smaller than 0 in y
[tempInd, tempInd2] = find(subpixMaxima(:,2,:) < 0);
subpixMaxima(tempInd,1:3,tempInd2) = 0;

%view pixel resolution maxima
figure
scatter3(rM,cM,sM,'.')

%view subpixel resolution maxima
figure
for i = 1:noImgs
    if min(find(sM == i)) > 0
        scatter3(subpixMaxima(:,1,i),subpixMaxima(:,2,i),subpixMaxima(:,3,i),'.')
        hold on
    end
end

%% Linking objects to pillars
% The plan is to create several metrics for determining whether a local 2D
% maxima belongs to a 'pillar' group of maxima by comparing the xy distance
% between the object of interest and the nearest neighbors on frames before
% and after.
close all

% Set a maximum linking distance in microns that any object can still be
% considered part of a pillar. Smaller values will speed up code.
maxLinkDistance = 1;
maxLD = maxLinkDistance/pixelSize;

% Set a maximum number of frames to look for a linked object before giving
% up (maxJumpDistance)
maxJD = 4;

% Find number of objects per frame in subpixMaxima
for i = 1:size(subpixMaxima,3)
    [lastNZElement,~] = find(subpixMaxima(:,1,i),1,'last');
    twoDimMaxIndSub(i,1) = lastNZElement;
end

%closest neighbor in frame above
noPillars = 0;
ignoreM = 1;
noProblemPillars = 0;
for i = 1:size(subpixMaxima,3)
    for j = 1:(twoDimMaxIndSub(i,1))
        
        %if on first frame, assign the current object a new pillar index
        if i == 1
            noPillars = noPillars +1;
            subpixMaxima(j,6,i) = noPillars;
        end
        ignoreN = 0; %end the loop early if this value changes
        for n = 1:maxJD %check all frames within jump distance range
            if ignoreN ==0 && i+n <= size(subpixMaxima,3)
                clear tempDistances
                
                %calulate distances from maxima in frame above to current object
                tempDistances = sqrt((subpixMaxima(1:twoDimMaxIndSub(i,1),1,i+n)-subpixMaxima(j,1,i)).^2 +(subpixMaxima(1:twoDimMaxIndSub(i,1),2,i+n)-subpixMaxima(j,2,i)).^2);
                
                %name the object in the next frame with the least distance from
                %current object its "nearest neighbor"
                [nearUpNeighbor,~] = find(tempDistances==min(tempDistances));
                
                %if the nearest neighbor falls within the max linking distance,
                if tempDistances(min(nearUpNeighbor),1) < maxLD
                    
                    %store the nearest up neighbor as as the fourth index in
                    %subpixMaxima
                    subpixMaxima(j,4,i) = min(nearUpNeighbor);
                    
                    %store the distance to the nearest up neighbor as index 5
                    subpixMaxima(j,5,i) = tempDistances(min(nearUpNeighbor),1);
                    subpixMaxima(j,8,i) = n;
                    ignoreN = 1;
                    
                    %assign a later pillar index if the current object has been
                    %previously linked
                    if subpixMaxima(j,6,i) > 0
                        subpixMaxima(subpixMaxima(j,4,i),6,i+n) = subpixMaxima(j,6,i);
                        
                        %if no match was found, assign a new pillar
                    else
                        noPillars = noPillars +1;
                        subpixMaxima(j,6,i) = noPillars;
                        subpixMaxima(subpixMaxima(j,4,i),6,i+n) = subpixMaxima(j,6,i);
                    end
                end
            end
            
            %%%%Try to assign pillar index to current pillar and pillar ahead!
            
            %%%%%%%%%%%%%%%%For Finding Pillars assigned in Previous Frames(no longer needed)
            %             if i>1
            %                 ignoreM = 0;
            %                 for m = 1:maxJD % looking back for previous links
            %                     if m<i && ismember(j,subpixMaxima(:,4,i-m))==1 && ignoreM==0
            %                         clear tempInd
            %                         tempInd = find(subpixMaxima(:,4,i-m)==j);
            %                         for tempCount = 1:size(tempInd,1)
            %                             ignoreTempCount = 0;
            %                             if subpixMaxima(tempInd(tempCount,1),8,i-m) == m && ignoreTempCount ==0
            %                                 subpixMaxima(j,6,i) = subpixMaxima(find(subpixMaxima(:,4,i-m)==j,1,'first'),6,i-m);
            %                                 subpixMaxima(j,7,i) = size(find(subpixMaxima(:,4,i-m)==j),1);
            %                                 ignoreM = 1;
            %                                 ignoreTempCount = 1;
            %                                 if m>1
            %                                     noProblemPillars = noProblemPillars+1
            %                                     problemPillars(noProblemPillars,1) = subpixMaxima(j,6,i);
            %                                 end
            %                             end
            %                         end
            %                     end
            %                 end
            %             end
            %%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

%% WIP second pass for correcting pillar errors
% l=0;
% searchWindowRadius = 1; %microns
% sWR = searchWindowRadius/pixelSize; %pixels
% for i = maxJD:size(subpixMaxima,3)-1
%     for j = 1:(twoDimMaxIndSub(i+1,1))
%         if size(find(subpixMaxima(:,4,i)==j),1)>1
%             culprits = find(subpixMaxima(:,4,i)==j);
%             for k = 1:size(culprits,1)
%                 if size(find(subpixMaxima(:,4,i-1)==culprits(k,1)),1)==0%if there is a hole before current trajectory
%                 %try to fill that hole and fix the pillar
%                 l = l+1
%                 subpixMaxima(culprits(k,1),6,i)
%                 clear Xs Ys
%                 sWX(1,1) = subpixMaxima(culprits(k,1),1,i)-sWR;
%                 sWX(1,2) = subpixMaxima(culprits(k,1),1,i)+sWR;
%                 sWY(1,1) = subpixMaxima(culprits(k,1),2,i)-sWR;
%                 sWY(1,2) = subpixMaxima(culprits(k,1),2,i)+sWR;
%                 [Xs(:,1),Xs(:,2)] = find(subpixMaxima(:,1,i-maxJD:i)>sWX(1,1) & subpixMaxima(:,1,i-maxJD:i)<sWX(1,2));
%                 [Ys(:,1),Ys(:,2)] = find(subpixMaxima(:,2,i-maxJD:i)>sWY(1,1) & subpixMaxima(:,2,i-maxJD:i)<sWY(1,2));
%                 for m = 1:size(Xs,1)
%                     tempXs = Xs(m,1:2);
%
%                 end
%
%                 else %if there is not a hole
%                 associates = find(subpixMaxima(:,4,i-1)==culprits(k,1));
%                 v1(1,1) = subpixMaxima(culprits(k,1),1,i) - subpixMaxima(associates(1,1),1,i-1);
%                 v1(1,2) = subpixMaxima(culprits(k,1),2,i) - subpixMaxima(associates(1,1),2,i-1);
%                 v2(1,1) = subpixMaxima(j,1,i+1) - subpixMaxima(culprits(k,1),1,i);
%                 v2(1,2) = subpixMaxima(j,2,i+1) - subpixMaxima(culprits(k,1),2,i);
%                 culprits(k,2) = sqrt((v2(1,1)-v1(1,1))^2 + (v2(1,2)-v1(1,2))^2);
%                 end
%             end
%         end
%     end
% end
%%
clear pillarBook
pillarBook = zeros(noImgs,5,noPillars);
for i = 1:noPillars
    [row,frame] = find(subpixMaxima(:,6,:)==i);
    for j = 1:size(row,1)
        pillarBook(j,1,i) = subpixMaxima(row(j,1),1,frame(j,1));
        pillarBook(j,2,i) = subpixMaxima(row(j,1),2,frame(j,1));
        pillarBook(j,3,i) = subpixMaxima(row(j,1),3,frame(j,1));
        pillarBook(j,4,i) = i;
        if round(pillarBook(j,1,i)) > 0 && round(pillarBook(j,2,i))>0 && round(pillarBook(j,3,i))>0
            pillarBook(j,5,i) = roiImgs(round(pillarBook(j,1,i)),round(pillarBook(j,2,i)),pillarBook(j,3,i));
        else
            pillarBook(j,5,i) = 0;
        end
    end
end
createExcelForTrajectories(pillarBook);
%% 3D Scatterplot of points color coded by pillar
figure
for j = 1:noPillars
    scatter3(pillarBook(:,1,j),pillarBook(:,2,j),pillarBook(:,3,j),'.')
    hold on
end
hold off
%% 3D Plot of points color coded by pillar and connected
figure
for j = 1:noPillars
    clear tempPillar
    first = find(pillarBook(:,1,j),1,'first');
    last = find(pillarBook(:,1,j),1,'last');
    tempPillar = pillarBook(first:last,:,j);
    plot3(tempPillar(:,1),tempPillar(:,2),tempPillar(:,3))
    hold on
end
hold off
%% 3D Plot of pillars thought to have issues
figure
for j = 1:noProblemPillars
    clear tempPillar
    first = find(pillarBook(:,1,problemPillars(j,1)),1,'first');
    last = find(pillarBook(:,1,problemPillars(j,1)),1,'last');
    tempPillar = pillarBook(first:last,:,problemPillars(j,1));
    plot3(tempPillar(:,1),tempPillar(:,2),tempPillar(:,3))
    hold on
end
hold off

%% Kovesi's function subpix3d
% % Obtain vectors with coordinates for x,y,z positions of local maxima with
% % pixel resolution
% [rM,cM,sM] = ind2sub(size(localMaxima3D),find(localMaxima3D == 1));
%
% % Find subpixel maxima based on initial guesses from imregionalmax on the
% % original images
% [rsM,csM,ssM] = subpix3d(rM,cM,sM,cropImages2);
%
% % 3D plot subpixel local maxima
% figure
% scatter3(rsM,csM,ssM,'.')

%% Dividing maxima by plane
% %Find the center of each z plane of dots based on histogram of local 3D
% %maxima
% [~,zCenters] = find(imregionalmax(imgaussfilt(histcounts(sM,noImgs),3)));
%
% % Find the spacing (# of frames) between each plane
% zSpacing = zeros(1,size(zCenters,2)-1);
% for i = 1:(size(zCenters,2)-1)
%     zSpacing(1,i) = zCenters(1,i+1)-zCenters(1,i);
% end
%
% % Choose a 'tail' size around plane centers for associating local maxima
% % above and below the plane of interest.
% zTails = round((mean(zSpacing(1,:)))/4);
%
% % Assign local maxima to planes
% zPlaneIndices = zeros(size(zCenters,2),2);
% for i = 1:size(zCenters,2)
%     zPlaneIndices(i,1) = min(find(sM>(zCenters(1,i)-zTails)));
%     zPlaneIndices(i,2) = max(find(sM<(zCenters(1,i)+zTails)));
%     zPlaneIndices(i,3) = zPlaneIndices(i,2)-zPlaneIndices(i,1);
% end
%
% zSortedMaxima = zeros(max(zPlaneIndices(:,3)),3,size(zCenters,2));
% for i = 1:size(zCenters,2)
%     zSortedMaxima(1:zPlaneIndices(i,3)+1,1,i) = rM(zPlaneIndices(i,1):zPlaneIndices(i,2),1);
%     zSortedMaxima(1:zPlaneIndices(i,3)+1,2,i) = cM(zPlaneIndices(i,1):zPlaneIndices(i,2),1);
%     zSortedMaxima(1:zPlaneIndices(i,3)+1,3,i) = sM(zPlaneIndices(i,1):zPlaneIndices(i,2),1);
% end
%
% zSortedMaximaFixed = zSortedMaxima;
% for i = 1:size(zCenters,2)
%     temp = zeros(zPlaneIndices(i,3),3);
%     temp(1:zPlaneIndices(i,3),:) = zSortedMaxima(1:zPlaneIndices(i,3),:,i);
%     fitobject{i} = fit([temp(:,1),temp(:,2)],temp(:,3),'poly11');
%     plot(fitobject{i},[temp(:,1),temp(:,2)],temp(:,3))
%     zSortedMaximaFixed(:,3,i) = fitobject{i}(zSortedMaximaFixed(:,1,i),zSortedMaximaFixed(:,1,i));
%     hold on
% end
% hold off
%
%% Predict New Local 3D maxima by Plane
%
% for i = 1:size(zCenters,2)
%     temp = zeros(zPlaneIndices(i,3),3);
%     temp(1:zPlaneIndices(i,3),:) = zSortedMaximaFixed(1:zPlaneIndices(i,3),:,i);
%     plot(fitobject{i},[temp(:,1),temp(:,2)],temp(:,3))
%     hold on
% end
% hold off
%
%%
% % Find subpixel maxima based on initial guesses from plane fits on the
% % original images
% zSortedMaximaFixedSubpix = zeros(size(zSortedMaximaFixed,1),3,size(zCenters,2));
% for i = 1:size(zCenters,2)
% [rsM2,csM2,ssM2] = subpix3d(zSortedMaximaFixed(:,1,i),zSortedMaximaFixed(:,2,i),round(zSortedMaximaFixed(:,3,i)),cropImages2);
% zSortedMaximaFixedSubpix(1:size(rsM2,2),1,i) = rsM2;
% zSortedMaximaFixedSubpix(1:size(rsM2,2),2,i) = csM2;
% zSortedMaximaFixedSubpix(1:size(rsM2,2),3,i) = ssM2;
% end
%
% for i = 1:size(zCenters,2)
%     plot(fitobject{i},[zSortedMaximaFixedSubpix(:,1,i),zSortedMaximaFixedSubpix(:,2,i)],zSortedMaximaFixedSubpix(:,3,i))
%     hold on
% end
% hold off
% figure
% for i = 1:size(zCenters,2)
% scatter3(zSortedMaximaFixedSubpix(:,1,i),zSortedMaximaFixedSubpix(:,2,i),zSortedMaximaFixedSubpix(:,3,i),'.')
% hold on
% end
% hold off