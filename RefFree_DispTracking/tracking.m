function tracking(directory)
if nargin == 1
    cd(directory);
    auto = 1;
else
    auto =0;
end
set(0,'defaultfigurecolor',[1 1 1])
% Ensures the path to necessary functions is available to the rest of the
% script
addpath(genpath('Tracking Functions'));
addpath(genpath('Kovesi Filters'));

%% Get Images and Metadata
disp('Creating Pre-Processed Images.')
if auto == 0
experiment = Experiment;
elseif auto == 1
[~,stackfile] = loadSpecific('Original.tif','*.tif','Did not find your file! Load a stack!')    
[~,Transfile] = loadSpecific('PlaceHolder.tif','*.tif','Did not find your file! Load a Transmitted Image!');
[~,Fluorfile] = loadSpecific('PlaceHolder.tif','*.tif','Did not find your file! Load a Transmitted Image!');
    
experiment = Experiment(stackfile,Transfile,Fluorfile);
end
% Scaling is the same in X and Y; convert from meters to microns
pixelSize = experiment.metadata.scalingX*1000000;
scaleFactor = (experiment.metadata.scalingX*1000000)/0.1625;
disp('Done.')
disp(scaleFactor)
disp(pixelSize)

%% Initial Image cropping (cannot be undone except by restarting script)
if auto == 0
[experiment.images,experiment.fluorImg,experiment.cellImg,roiBounds1] = experiment.cropImgs;
end

%% Final Pre-processing Before Finding Local Maxima
clear roiCell;
[experiment.ppOptions,experiment.masks] = experiment.preprocess(auto);
[roiImgs,roiMasks,roiCell,roiBounds,roiZeros,redoCheck] = experiment.cropImgs2(auto);

while redoCheck == 1
    experiment = Experiment(experiment.images,experiment.cellImg,experiment.fluorImg,experiment.metadata);
    [experiment.ppOptions,experiment.masks] = experiment.preprocess;
    pixelSize = experiment.metadata.scalingX*1000000;
    scaleFactor = (experiment.metadata.scalingX*1000000)/0.165;
    disp('Scale Factor:')
    disp(scaleFactor)
    [roiImgs,roiMasks,roiCell,roiBounds,roiZeros,redoCheck] = experiment.cropImgs2;
end
%% Remove Black Regions
%This should help with image stacks that are not completely filled with dot
%array. Without this, the object detection protocol detects false objects
%in the noise of an unpatterned region.
% if ismember(5,experiment.ppOptions{1}) ==0 && ismember(6,experiment.ppOptions{1}) ==0
% close all
% se = strel('disk',5);
% se2 = strel('disk',10);
% empty = roiMasks(:,:,1)>0;
% emptyMask = imerode(bwareaopen(imclose(imdilate(empty,se),se2),1000),se);
% emptyRegion = figure;
% imshow(emptyMask,[])
% filePath = cd;
% savefile = [filePath '\Mask of Empty Regions.tif'];
% export_fig(emptyRegion,savefile,'-native');
% else
emptyMask = ones([size(roiMasks(:,:,1),1) size(roiMasks(:,:,1),2)]);
% end

%% Viewing Pillars as a Frame Weighted Z-Projection
%Frames are thresholded and projected through-Z. Pixels appearing in later
%frames appear brightest.
for i = 1:size(roiMasks,3)
    temp = roiMasks(:,:,i).*emptyMask;
    low = mean(mean(temp(temp>0)));
    
    high = max(max(roiMasks(:,:,i)));
    temp(temp>0) = high*(1);
    temp = (temp/high)*(i^4/size(roiMasks,3)^4);
    
    roiMasks2(:,:,i) = temp;
end
maxMasks = permute(max(permute(roiMasks2,[3 1 2])),[2 3 1]);
pillarView = figure;
imshow(maxMasks,[])
filePath = cd;
savefile = [filePath '\Tracking_pillarView.tif'];
export_fig(pillarView,savefile,'-native');

%% Sub-pixel Object Detection Using Kilfoil 'feature2D' Function
clear subpixMaxima3
%[r] = feature2D(img,lambda,w,masscut,Imin)
for i = 1:size(roiImgs,3)
    clear temp
    currentImg = roiMasks(:,:,i).*emptyMask;
    temp = feature2D(currentImg,1,2,2000,50);
    if i == 1
        temp(:,6) = i;
        subpixMaxima3 = temp;
    else
        temp(:,6) = i;
        subpixMaxima3 = cat(1,subpixMaxima3,temp);
    end
end

%% Object Detection Plot
detections = figure;
imshow(roiZeros)
hold on
map = brewermap(max(subpixMaxima3(:,6)),'*Spectral');

for i = 1:max(subpixMaxima3(:,6)) %size(subpixMaxima,3)-10
    [tempFrame,~] = find(subpixMaxima3(:,6) == i);
    scatter3(subpixMaxima3(tempFrame,1),subpixMaxima3(tempFrame,2),subpixMaxima3(tempFrame,6),'.','SizeData',100,'MarkerFaceColor',map(i,1:3),'MarkerEdgeColor',map(i,1:3))
end
hold off
filePath = cd;
savefile = [filePath '\Tracking_Unlinked Detections.tif'];
export_fig(detections,savefile,'-native');
%%
% Linking objects to pillars (i.e. Through the z dim)
% The plan is to create several metrics for determining whether a local 2D
% maxima belongs to a 'pillar' group of maxima by comparing the xy distance
% between the object of interest and the nearest neighbors on frames before
% and after.
close all
disp('Linking Dots Between Frames.')
% Set a maximum linking distance in microns that any object can still be
% considered part of a pillar. Smaller values will speed up code.
maxLinkDistance = 1;
maxLD = maxLinkDistance/pixelSize;
disp(['Max Link Distance (Microns): ',num2str(maxLinkDistance)])
% Set a maximum number of frames to look for a linked object before giving
% up (maxJumpDistance)
maxJD = 15;
disp(['Max Jump Distance (Frames): ',num2str(maxJD)])

%% Linking using trackmem function (kilfoil code)
disp('Linking Pillars')
[lub] = trackmem(subpixMaxima3,maxLD,2,1,maxJD);
disp('Done Linking Pillars')
% Cleaning tracks and fixing incorrect matches
maxAi = 45;
maxAc = 90;
lub(:,8) = 1:1:size(lub,1);
maxD = 5/pixelSize;

%This sequence preps "lub" for the refreshData fxn which calculates useful
%properties of individual pillars, including a quality score to determine
%whether a particular pillar has an issue.
lub = sortrows(lub,[7 6]);
numPillars = max(lub(:,7));
[lub] = refreshData(lub,maxAi,maxAc,numPillars,maxMasks);
clear problems
problems = findProbs(lub);

%View pillars and overlay "Problems"
% plotTracks(lub,maxMasks)
% hold on
% scatter3(lub(problems,1),lub(problems,2),lub(problems,6))

% Automated Tracking Corrections
%Find tracks that are shorter than they should be and use their history to
%fix them.
nSize = 20; % max neighbor distance
tic
disp(['Attempting to Autofix Tracking Errors'])

for j = 1 % number of iterations
    clear pts
    problems(:,3) = 0; % third column will track if problem object has been addressed
    for L = 1:size(problems,1) % For each problem pillar
        replace = 0;
        %If the problem has not been addressed:
        if problems(L,3) ~= 1
            clear cntrPt pPilIdx
            pPilIdx = problems(L,1);
            pPilIdcs = find(lub(:,7)==lub(problems(L,1),7) & lub(:,6)>=lub(problems(L,1),6)); %Find all objects with common pillar with current object and above current object
            pPilIdcs2 = find(lub(problems(:,1),7)==lub(problems(L,1),7)); % Also note ALL members of current objects pillar that are also problems
            cntrPt = round(lub(pPilIdx,1:2)); %Find an XY coordinate describing the *center* of current objects pillar
            
            clear tempInd1
            tempInd1 = find(lub(problems(pPilIdcs2),7)==lub(pPilIdx,7));
            if size(tempInd1,1)>1
                tempPillar = lub(problems(pPilIdcs2),:);
                tempPillar = sortrows(tempPillar,6);
                useTI = 0;
            else
                tempInd1 = find(lub(:,7)==lub(pPilIdx,7) & lub(:,6)>=lub(pPilIdx,6));
                tempPillar = lub(tempInd1,:);
                tempPillar = sortrows(tempPillar,6);
                useTI = 1;
            end
            
            clear potMatch
            potMatch = 0;
            if size(tempPillar,1) > 1
                nghbrs = unique(lub((lub(:,1)>(cntrPt(1,1)-nSize)&lub(:,1)<(cntrPt(1,1)+nSize)&lub(:,2)>(cntrPt(1,2)-nSize)&lub(:,2)<(cntrPt(1,2)+nSize)),7));
                %---Approximate Current Pillars trajectory
                %Create Trajectory vector of length 1um
                length = 1; %microns
                pixel = experiment.metadata.scalingX*1000000; %pixel size in microns
                
                if tempPillar(end,12) > 1 && tempPillar(end,27)>3
                    [ppd,perp] = meanTraj(tempPillar,length,pixel);
                else
                    clear paraM perpM
                    for i = 1:size(nghbrs,1)
                        clear tempPillar
                        tempPillar = lub(lub(:,7)==nghbrs(i,1),:);
                        [para,perp] = meanTraj(tempPillar,length,pixel);
                        if i == 1
                            paraM = para;
                            perpM = perp;
                        else
                            paraM = cat(1,paraM,para);
                            perpM = cat(1,perpM,perp);
                        end
                    end
                    ppd = mean(paraM,1);
                    rowScale = sqrt(((length/pixel)^2)/(ppd(1,1)^2+ppd(1,2)^2));
                    ppd = ppd *rowScale;
                    perp = mean(perpM,1);
                    rowScale = sqrt(((length/pixel)^2)/(perp(1,1)^2+perp(1,2)^2));
                    perp = perp *rowScale;
                end
                
                %Create Rectangular Search window to find potential matches
                pts(1,1:2,L) = lub(pPilIdx,1:2) - 1*perp(1,1:2) + 0*ppd(1,1:2);
                pts(2,1:2,L) = lub(pPilIdx,1:2) + 1*perp(1,1:2) + 0*ppd(1,1:2);
                pts(3,1:2,L) = lub(pPilIdx,1:2) + 1*perp(1,1:2) + 3*ppd(1,1:2);
                pts(4,1:2,L) = lub(pPilIdx,1:2) - 1*perp(1,1:2) + 3*ppd(1,1:2);
                
                
                %---Use trajectory information to determine best match from
                %neighbor pillars
                
                %Index neighbors and neighbor pillars
                
                
                
                for i = 1:size(nghbrs,1)
                    if size(find(lub(:,7)==nghbrs(i,1) & lub(:,6)<lub(pPilIdx,6) & lub(:,30)~=1,1,'last'),1)>0
                        nghbrs(i,2) = find(lub(:,7)==nghbrs(i,1) & lub(:,6)<lub(pPilIdx,6) & lub(:,30)~=1,1,'last');
                        if inpolygon(lub(nghbrs(i,2),1),lub(nghbrs(i,2),2),pts(1:4,1,L),pts(1:4,2,L))
                            if  max(lub(lub(:,7)==nghbrs(i,1),6)) < lub(pPilIdx,6) % The case where a solution exists where the uppermost object in solution is lower than lowermost object in pillar
                                if potMatch == 0
                                    potMatch = nghbrs(i,2);
                                    replace = 0;
                                else
                                    potMatch = cat(1,potMatch,nghbrs(i,2));
                                    replace = 0;
                                end
                            elseif max(lub(lub(:,7)==nghbrs(i,1),6)) >= lub(pPilIdx,6) && useTI~=1 % The case where a solution exists but the solution currently includes its own errors
                                if size(pPilIdcs2,1) >1
                                    if potMatch == 0
                                        potMatch = nghbrs(i,2);
                                        replace = 1;
                                    else
                                        potMatch = cat(1,potMatch,nghbrs(i,2));
                                        replace = 1;
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
            
            %DEBUGGING!!!!!
            %             if potMatch(1,1)~=0 %&& j>1
            %                 w = figure;
            %                 hold on
            %                 if size(potMatch,1)>1
            %                     diffs = (lub(pPilIdx,1:2)-lub(potMatch,1:2)).^2;
            %                     diffs2 = sqrt(diffs(:,1)+diffs(:,2));
            %                     closeIdx = find(diffs2 == min(diffs2),1,'first');
            %                 else
            %                     closeIdx = 1;
            %
            %                 end
            %                 axis([cntrPt(1,1)-fSizeXmin cntrPt(1,1)+fSizeXmax cntrPt(1,2)-fSizeYmin cntrPt(1,2)+fSizeYmax])
            %                 plot3(tempPillar(:,1),tempPillar(:,2),tempPillar(:,6))
            %                 plot(pts(:,1),pts(:,2))
            %                 scatter3(lub(potMatch(closeIdx),1),lub(potMatch(closeIdx),2),lub(potMatch(closeIdx),6))
            %                 tarPil = find(lub(:,7)== lub(potMatch(closeIdx),7));
            %                 plot3(lub(tarPil,1),lub(tarPil,2),lub(tarPil,6))
            %                 waitfor(w)
            %             end
            %DEBUGGING!!!!! END
            
            
            
            if size(potMatch,1)>1
                for p = 1:size(potMatch,1)
                    if lub(potMatch,7) == lub(pPilIdx,7)
                        potMatch(p,:) =[];
                    end
                end
            end
            if size(potMatch,1)>1
                diffs = (lub(pPilIdx,1:2)-lub(potMatch,1:2)).^2;
                diffs2 = sqrt(diffs(:,1)+diffs(:,2));
                closeIdx = find(diffs2 == min(diffs2),1,'first');
                if replace == 1
                    replaceIdc = find(lub(:,7)==lub(potMatch(closeIdx),7) & lub(:,6) >= lub(pPilIdx,6));
                    lub(replaceIdc,7) = max(lub(:,7)+1);
                end
                if useTI == 0
                    lub(problems(pPilIdcs2(:)),7) = lub(potMatch(closeIdx),7);
                elseif useTI == 1
                    lub(tempInd1,7) = lub(potMatch(closeIdx),7);
                end
                problems(pPilIdcs2,3) = 1;
            elseif potMatch~=0
                if replace == 1
                    replaceIdc = find(lub(:,7)==lub(potMatch(1),7) & lub(:,6) >= lub(pPilIdx,6));
                    lub(replaceIdc,7) = max(lub(:,7)+1);
                end
                if useTI == 0
                    lub(problems(pPilIdcs2(:)),7) = lub(potMatch(1,1),7);
                elseif useTI == 1
                    lub(tempInd1,7) = lub(potMatch(1,1),7);
                    
                end
                problems(pPilIdcs2,3) = 1;
            end
            
        end
    end
    
    lub = sortrows(lub,[7 6]);
    [lub] = refreshData(lub,maxAi,maxAc,numPillars,maxMasks);
    problems = findProbs(lub);
end
disp(['Done Attempting to Fix Errors at '  num2str(toc) ' seconds'])


% Automated Problem Pillar Selection with Manual Track Adjustment
%This code will allow the user to select a region from an image to address
%the "problems" in that region. Many errors occur because the input
%parameters for the tracking function are pretty conservative, but can be
%easily fixed manually.
if auto == 0
w = questdlg('Would you like to correct tracking errors?',...
    'Manual Adjustments (Optional)','Yes','No','Yes');
waitfor(w);
close()
if strcmp(w,'Yes')==1
    %Allow the user to first choose an area to be fixed
    fixProblems = 1;
    while fixProblems == 1
        plotTracks(lub,ones(size(maxMasks)))
        hold on
        %pOffAngle = find(problems(:,2));
        pBackward = find(problems(:,3));
        pEndAngle = find(problems(:,4));
        pSize = find(problems(:,5));
        
        %scatter3(lub(problems(pOffAngle,1),1),lub(problems(pOffAngle,1),2),lub(problems(pOffAngle,1),6),'black')
        scatter3(lub(problems(pBackward,1),1),lub(problems(pBackward,1),2),lub(problems(pBackward,1),6),'r')
        scatter3(lub(problems(pEndAngle,1),1),lub(problems(pEndAngle,1),2),lub(problems(pEndAngle,1),6),'g')
        scatter3(lub(problems(pSize,1),1),lub(problems(pSize,1),2),lub(problems(pSize,1),6),'b')
        
        for i = 1:size(pts,3)
            plot(pts(:,1,i),pts(:,2,i))
        end
        
        
        w = questdlg('Select a region containing incorrect tracking. Double-Click to confirm selection.',...
            'Select problem area','Okay','Okay');
        [~,pBds]=imcrop;
        close()
        pBds = round(pBds);
        pBds(1,3) = pBds(1,3)+pBds(1,1);
        pBds(1,4) = pBds(1,4)+pBds(1,2);
        
        for j = 1
            problems(:,3) = 0;
            for l = 1:size(problems,1)
                if lub(problems(l,1),1)>pBds(1,1) && lub(problems(l,1),1)<pBds(1,3) && lub(problems(l,1),2)>pBds(1,2) && lub(problems(l,1),2)<pBds(1,4)
                    if problems(l,3) ~= 1
                        clear cntrPt pPilIdx
                        pPilIdx = problems(l,1);
                        pPilIdcs = find(lub(:,7)==lub(problems(l,1),7) & lub(:,6)>=lub(problems(l,1),6));
                        pPilIdcs2 = find(lub(problems(:,1),7)==lub(problems(l,1),7));
                        cntrPt = round(lub(pPilIdx,1:2));
                        fSize = 50; % Distance from center point to make frame
                        nSize = 20; % max neighbor distance
                        
                        %determine if object is near the edge of the main frame
                        if size(maxMasks,2)-cntrPt(1,1) < fSize
                            fSizeXmin = fSize;fSizeXmax = (size(maxMasks,2)-cntrPt(1,1))-1;
                        elseif cntrPt(1,1) <= fSize
                            fSizeXmin = cntrPt(1,1)-1;fSizeXmax = fSize;
                        else
                            fSizeXmin = fSize;fSizeXmax = fSize;
                        end
                        
                        if size(maxMasks,1)-cntrPt(1,2) < fSize
                            fSizeYmin = fSize;fSizeYmax = (size(maxMasks,1)-cntrPt(1,2))-1;
                        elseif cntrPt(1,2) <= fSize
                            fSizeYmin = cntrPt(1,2)-1;fSizeYmax = fSize;
                        else
                            fSizeYmin = fSize;fSizeYmax = fSize;
                        end
                        
                        pPF = figure('Position',[100 100 800 800]);
                        handles.axes1 = axes( ...
                            'Units','normalized', ...
                            'Position',[0 0 1 1]);
                        imshow(maxMasks(cntrPt(1,2)-fSizeYmin:cntrPt(1,2)+fSizeYmax,cntrPt(1,1)-fSizeXmin:cntrPt(1,1)+fSizeXmax))
                        hold on
                        %pbd = uicontrol('style','pushbutton','units','pixels','position',[120,5,70,20],'string','delete','Callback',@pbDelete);
                        pbs = uicontrol('style','pushbutton','units','pixels','position',[120,25,70,20],'string','skip','Callback',@pbSkip);
                        
                        
                        %Plot Problem
                        scatter3(fSizeXmin,fSizeYmin,lub(pPilIdx,6),'g')
                        
                        %Plot neighbors and neighbor pillars
                        nghbrs = unique(lub((lub(:,1)>(cntrPt(1,1)-nSize)&lub(:,1)<(cntrPt(1,1)+nSize)&lub(:,2)>(cntrPt(1,2)-nSize)&lub(:,2)<(cntrPt(1,2)+nSize)),7));
                        for i = 1:size(nghbrs,1)
                            if size(find(lub(:,7)==nghbrs(i,1) & lub(:,6)<lub(pPilIdx,6),1,'last'),1)>0
                                nghbrs(i,2) = find(lub(:,7)==nghbrs(i,1) & lub(:,6)<lub(pPilIdx,6),1,'last');
                                scatter3(lub(nghbrs(i,2),1)-(cntrPt(1,1)-fSizeXmin),lub(nghbrs(i,2),2)-(cntrPt(1,2)-fSizeYmin),lub(nghbrs(i,2),6),'r')
                                trckNum = num2str(nghbrs(i,1));
                                tempInd1 = find(lub(:,7)==nghbrs(i,1));
                                if size(tempInd1,1)>0
                                    tempPillar = lub(tempInd1,:);
                                    tempPillar = sortrows(tempPillar,6);
                                    plot3(tempPillar(:,1)-(cntrPt(1,1)-fSizeXmin),tempPillar(:,2)-(cntrPt(1,2)-fSizeYmin),tempPillar(:,6))
                                    %plot3(tempPillar(:,24),tempPillar(:,25),tempPillar(:,6))
                                end
                                pb(i) = uicontrol('style','pushbutton','units','pixels','position',[40,5+(20*i),70,20],'string',num2str(trckNum),'ButtonDownFcn',@pbPillarChoice,'Callback',@pbPillarChoice);
                                
                                trckText = strcat('\leftarrow ',trckNum);
                                text(lub(nghbrs(i,2),1)-(cntrPt(1,1)-fSizeXmin),lub(nghbrs(i,2),2)-(cntrPt(1,2)-fSizeYmin),lub(nghbrs(i,2),6),trckText,'Color','red')
                            else
                                nghbrs(i,2) = 0;
                            end
                        end
                        %Plot Problem's Current Pillar
                        tempInd1 = find(lub(:,7)==lub(pPilIdx,7));
                        if size(tempInd1,1)>0
                            tempPillar = lub(tempInd1,:);
                            tempPillar = sortrows(tempPillar,6);
                            plot3(tempPillar(:,1)-(cntrPt(1,1)-fSizeXmin),tempPillar(:,2)-(cntrPt(1,2)-fSizeYmin),tempPillar(:,6),'g')
                        end
                        pPF.UserData = zeros(3,1);
                        waitfor(pPF,'UserData')
                        choice = pPF.UserData;
                        if choice(1,1) == -1
                            problems(l,4) = 1;
                            problems(pPilIdcs2,4) = 1;
                        elseif choice(1,1) == -2
                        else
                            %         if size(find(lub(:,7)==choice & lub(:,6)>=lub(pPilIdx,6)),1)>0
                            %         replace = find(lub(:,7)==choice & lub(:,6)>=lub(pPilIdx,6));
                            %         lub(replace,7) = max(lub(:,7)+1);
                            %         end
                            if choice(2,1) ~=0
                                disp('Right Click!')
                                find(lub(lub(:,7)==lub(choice(2,1)) & lub(:,6)>=lub(pPilIdx,6),7))
                                lub(lub(:,7)==choice(2,1) & lub(:,6)>=lub(pPilIdx,6),7) = max(lub(:,7))+1;
                                lub(pPilIdcs,7) = choice(2,1);
                            else
                                lub(pPilIdcs,7) = choice(1,1);
                            end
                            problems(pPilIdcs2,3) = 1;
                        end
                        close()
                    end
                end
            end
            
            %delete = find(problems(:,4)==1);
            %lub(problems(delete,1),:) = [];
            
            lub = sortrows(lub,[7 6]);
            [lub] = refreshData(lub,maxAi,maxAc,numPillars,maxMasks);
            problems = findProbs(lub);
        end
        w = questdlg('Fix another error?',...
            'Manual Adjustments (Optional)','Yes','No','Yes');
        if strcmp(w,'Yes') ~= 1
            fixProblems = 0;
        end
    end
end
end
%%
%Duplicate Point Correction
%Rationale: Some objects appear in two frames at once, typically at the
%junction between one ellipsoid and the next when there is large shear
%deformation. The following code should give the user the option to merge
%duplicate points (average their location) or choose a single one for a
%problem pillar.

%Simple Method
for i = 1:size(lub,1)-1
    if lub(i,6) == lub(i+1,6) && lub(i,7) == lub(i+1,7)
        lub(i,:) = (lub(i,:) + lub(i+1,:))./2;
        lub(i+1,:) = lub(i,:);
    end
end
[lub2,a,b] = unique(lub,'rows');
lub = lub2(b,:);

%% Remove Pillars Smaller than Threshold
sT = 10;
numPillars = unique(lub(:,7));
for i = 1:size(numPillars)
    tempInd1 = find(lub(:,7)==numPillars(i,1));
    if size(tempInd1,1) <= sT
        lub(tempInd1,:)=[];
    end
end

%% Creating Trajectories Input text file
disp('Creating Text File for trajectories.m.')
clear subpixMaxima2
%Organize subpixMaxima2 to be compatible with trajectories script
subpixMaxima2(:,1)=1:1:size(lub,1);
subpixMaxima2(:,2)=lub(:,1);
subpixMaxima2(:,3)=lub(:,2);
subpixMaxima2(:,4)=lub(:,6);
subpixMaxima2(:,5)=lub(:,7);
subpixMaxima2(:,6)=lub(:,3);
subpixMaxima2(:,7)=lub(:,4);
subpixMaxima2(:,8)=lub(:,5);

%Create the text file
txtFileName = 'trajectoriesInput.txt';

%Replace existing version
if exist(txtFileName,'file')
    delete(txtFileName)
end
dlmwrite(txtFileName,subpixMaxima2);


%% 2D Plot of points color coded by pillar and connected
detections2 = figure;
imshow(maxMasks)
hold on
clear tempInd1 tempInd2
for j = 1:max(lub(:,7))
    clear tempPillar
    tempInd1 = find(lub(:,7)==j,1,'first');
    tempInd2 = find(lub(:,7)==j,1,'last');
    if((tempInd2-tempInd1)>0)
        tempPillar = lub(tempInd1:tempInd2,:);
        plot3(tempPillar(:,1),tempPillar(:,2),tempPillar(:,6))
    end
end
hold off
savefile = [filePath '\Tracking_Linked Detections on pillarView.tif'];
export_fig(detections2,savefile,'-native');

%% Saving Parameters of Pre-Processing
%Only relevant if not using Kilfoil pre-processing
disp('Saving Pre-processing Parameters.')
parametersObj{1} = experiment.ppOptions;
parametersObj{2} = maxLinkDistance;
parametersObj{3} = maxJD;
parametersObj{4} = pixelSize;
parametersObj{5} = scaleFactor;
paraTxt = fopen('parameters.txt','wt');
fprintf(paraTxt,strcat('Original File Name: ', experiment.metadata.filename, '\n'));
p1Format = 'Remove Large?';
if ismember(2,parametersObj{1,1}{1,1})==1
    fprintf(paraTxt,p1Format);
    fprintf(paraTxt,' yes \n');
    p1Format = 'Remove Large Size(microns squared): %d \n';
    fprintf(paraTxt,p1Format,parametersObj{1,1}{1,4});
else
    fprintf(paraTxt,p1Format);
    fprintf(paraTxt,' no \n');
end

p1Format = 'Approximate Feature Diameter XY: %0.2f \n';
fprintf(paraTxt,p1Format,parametersObj{1,1}{1,2});

p1Format = 'Approximate Feature Diameter Z: %0.2f \n';
fprintf(paraTxt,p1Format,parametersObj{1,1}{1,3});

p1Format = 'Max Link Distance: %0.2f \n';
fprintf(paraTxt,p1Format,parametersObj{2});

p1Format = 'Max Jump Distance: %i \n';
fprintf(paraTxt,p1Format,parametersObj{3});

p1Format = 'Pixel Size: %f \n';
fprintf(paraTxt,p1Format,parametersObj{4});

p1Format = 'Scale Factor: %f \n';
fprintf(paraTxt,p1Format,parametersObj{5});

fclose(paraTxt);
%%
clear experiment roiCell roiImgs roiMasks roiMasks2
save tracking
disp('tracking.m is completed.')