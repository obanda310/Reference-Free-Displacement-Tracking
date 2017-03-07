function [SPM,numPillars] = LinkMaxima(SPM,maxLD,maxJD)
count = 0;
%SPM is subpixMaxima from tracking code
%TDMIS = 2-D Maximum Indices of SPM (number of objects on a page of SPM)

%SPM Structure
%(x,:,:) describes object ID
%(:,x,:) describes object parameter (the indices below)
%(:,:,x) describes object frame

%SPM Indices
%1 X
%2 Y
%3 Z (frame)
%4 index of nearest neighbor up within link distance
%5 XY distance to nearest neighbor above
%6 Pillar ID
%7 number of frames before current pillar that current pillar was
%previously linked to its pillar
%8 number of frames after current pillar to reach nearest neighbor
%9 Current objects place in pillar (to measure the number of members in a
%pillar)
%10 Original X
%11 Original Y
%12 Angle in degrees (ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));)
%13 Competitor's angle in degrees

% Find number of objects per frame in SPM
for i = 1:size(SPM,3)
    [lastNZElement,~] = find(SPM(:,1,i),1,'last');
    TDMIS(i,1) = lastNZElement;
end

%%%%%%%%%%%%%%%%%%%%%PILLAR LINKING START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up Outer Looping Sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numPillars = 0;
for i = 1:size(SPM,3)
    for j = 1:(TDMIS(i,1))
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Setting up First Pillars on First Frame
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %if on first frame, assign the current object a new pillar index
        %and continue
        if i == 1
            numPillars = numPillars +1; %increase pillar count
            SPM(j,6,i) = numPillars; %store pillar ID
            SPM(j,9,i) = 1; %record current pillar length (has to be 1 on first frame)
            SPM(j,10,i) = SPM(j,1,i); %store first frame of pillar XY
            SPM(j,11,i) = SPM(j,2,i); %redundant, but storing for consistency later
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Setting up Inner Looping Sequence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ignoreN = 0; %end the loop early if this value changes
        for n = 1:maxJD %check all frames within jump distance range
            if ignoreN ==0 && i+n <= size(SPM,3)
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Establish the preliminary best nearest neighbor
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                clear tempDistances nearUpNeighbors
                %calulate distances from maxima in frame above to current object
                tempDistances = sqrt((SPM(1:TDMIS(i,1),1,i+n)-SPM(j,1,i)).^2 +(SPM(1:TDMIS(i,1),2,i+n)-SPM(j,2,i)).^2);
                
                %name the object in the next frame with the least distance from
                %current object its "nearest neighbor"
                [nearUpNeighbors,~] = find(tempDistances<maxLD);
                
                
                skip = 0;
                if size(nearUpNeighbors,1) > 1 && SPM(j,9,i) > 2 %pillar length must be greater than 2 to have enough points to create an angle
                    %the case where there are multiple neighbors and
                    %pillar history exists
                    clear neighborAngles
                    for m = 1:size(nearUpNeighbors)
                        P0 = SPM(j,10:11,i);
                        P1 = SPM(j,1:2,i);
                        P2 = SPM(nearUpNeighbors(m,1),1:2,i+n);
                        neighborAngles(m,1) = (180/pi)*atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
                    end
                    clear neighborOptions
                    neighborOptions = find(neighborAngles==min(neighborAngles));
                    nearUpNeighbor = nearUpNeighbors(neighborOptions(1,1),1);
                    
                    disp = sqrt((SPM(j,1,i)-SPM(j,10,i)).^2 +(SPM(j,2,i)-SPM(j,11,i)).^2);
                    if min(neighborAngles) > 10 && disp>5 %(optional: ...&& tempDistances(nearUpNeighbor,1) > 2)
                        skip = 1;
                    end
                    
                    
                elseif size(nearUpNeighbors,1) == 1 && SPM(j,9,i) > 2
                    %the case where there is only a single neighbor and
                    %pillar history exists
                    nearUpNeighbor = nearUpNeighbors(1,1);
                    P0 = SPM(j,10:11,i);
                    P1 = SPM(j,1:2,i);
                    P2 = SPM(nearUpNeighbors(1,1),1:2,i+n);
                    neighborAngle = (180/pi)*atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0)); 
                    disp = sqrt((SPM(j,1,i)-SPM(j,10,i)).^2 +(SPM(j,2,i)-SPM(j,11,i)).^2);
                    if neighborAngle > 15 && disp>5 %(optional: ...&& tempDistances(nearUpNeighbor,1) > 2)
                        skip = 1;
                    end
                    
                else
                    %the simple case (not a lot of motion, no pillar history, so no angles)
                    tempDistances2 = tempDistances(nearUpNeighbors);
                    B = tempDistances2(:,1) == min(tempDistances2);
                    nearUpNeighbor = nearUpNeighbors(B,1);
                    if size(nearUpNeighbor,1)>1
                        clear nearUpNeighbor2
                        nearUpNeighbor2 = nearUpNeighbor(1,1);
                        clear nearUpNeighbor
                        nearUpNeighbor = nearUpNeighbor2;
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Pillar Assignment and Checks Against Previous Pillar Match Assignments
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %if the nearest neighbor falls within the max linking distance,
                if size(nearUpNeighbor,1)>0 && tempDistances(nearUpNeighbor)<maxLD && skip == 0 
                    
                    %store the nearest up neighbor as as the fourth index in
                    %SPM
                    SPM(j,4,i) = min(nearUpNeighbor);
                    
                    %store the distance to the nearest up neighbor as index 5
                    SPM(j,5,i) = tempDistances(min(nearUpNeighbor),1);
                    SPM(j,8,i) = n;
                    ignoreN = 1; %a match has been found the loop should proceed to next object after pillar assignment
                    
                    
                    
                    %if a pillar exists for current object and also it's
                    %nearest neighbor in following frame, determine
                    %whether current pillar is a better match. If it
                    %is, replace the nearest neighbor pillar
                    if SPM(j,6,i) > 0 && SPM(SPM(j,4,i),6,i+n) > 0
                        %Compare Angles first
                        
                        if SPM(SPM(j,4,i),9,i+n) > 2 && SPM(j,9,i) > 2
                            P0 = SPM(SPM(j,4,i),10:11,i+n);
                            P1 = SPM(j,1:2,i);
                            P2 = SPM(SPM(j,4,i),1:2,i+n);
                            SPM(SPM(j,4,i),13,i+n) = (180/pi)*atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
                            if SPM(SPM(j,4,i),13,i+n) < SPM(SPM(j,4,i),12,i+n)
                                SPM(SPM(j,4,i),6,i+n) = SPM(j,6,i);
                                SPM(SPM(j,4,i),7,i+n) = n;
                                SPM(SPM(j,4,i),9,i+n) = SPM(j,9,i)+1;
                            end
                        else
                            %Otherwise use pillar length
                            current = SPM(j,9,i);
                            [previousIdx,~] = find(SPM(:,6,(i+n)-SPM(SPM(j,4,i),7,i+n))==SPM(SPM(j,4,i),6,i+n));
                            previous = SPM(previousIdx(1,1),9,(i+n)-SPM(SPM(j,4,i),7,i+n));
                            if current > previous
                                SPM(SPM(j,4,i),6,i+n) = SPM(j,6,i);
                                SPM(SPM(j,4,i),7,i+n) = n;
                                SPM(SPM(j,4,i),9,i+n) = SPM(j,9,i)+1;
                            end
                        end
                        
                        %if a pillar exists for current object but not it's
                        %nearest neighbor in following frame, assign it's
                        %nearest neighbor the pillar
                    elseif SPM(j,6,i) > 0 && SPM(SPM(j,4,i),6,i+n) == 0
                        SPM(SPM(j,4,i),6,i+n) = SPM(j,6,i); %assign current pillar ID
                        SPM(SPM(j,4,i),7,i+n) = n; %store n to find nearest neighbor frame later if needed
                        SPM(SPM(j,4,i),9,i+n) = SPM(j,9,i)+1; %store pillar size as current size +1
                        SPM(SPM(j,4,i),10,i+n) = SPM(j,10,i);%store starting xy
                        SPM(SPM(j,4,i),11,i+n) = SPM(j,11,i);
                        
                        %if there are enough points to calculate an angle
                        %ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
                        if SPM(SPM(j,4,i),9,i+n) > 2
                            P0 = SPM(SPM(j,4,i),10:11,i+n);
                            P1 = SPM(j,1:2,i);
                            P2 = SPM(SPM(j,4,i),1:2,i+n);
                            SPM(SPM(j,4,i),12,i+n) = (180/pi)*atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
                        end
                        
                        
                        %if no pillar exists for current object or it's
                        %nearest neighbor in following frame, assign it
                        %a new pillar and assign it's nearest neighbor the
                        %pillar as well
                    elseif SPM(j,6,i) == 0 && SPM(SPM(j,4,i),6,i+n) == 0
                        numPillars = numPillars +1;
                        SPM(j,6,i) = numPillars;
                        SPM(SPM(j,4,i),6,i+n) = numPillars;
                        SPM(SPM(j,4,i),7,i+n) = n;
                        SPM(j,9,i) = 1; %current object is pillar size 1
                        SPM(SPM(j,4,i),9,i+n) = 2; %neighbor object makes pillar size 2
                        SPM(j,10,i) = SPM(j,1,i);%store starting xy
                        SPM(j,11,i) = SPM(j,2,i);
                        SPM(SPM(j,4,i),10,i+n) = SPM(j,10,i);%store starting xy
                        SPM(SPM(j,4,i),11,i+n) = SPM(j,11,i);
                        
%***** following doesn't happen under current logic*****                        
%                         %if no pillar exists for current object, but one
%                         %does exist for it's nearest neighbor in following
%                         %frame, create new pillar, and determine if current
%                         %pillar is a better match
%                         
%                     elseif SPM(j,6,i) == 0 && SPM(SPM(j,4,i),6,i+n) > 0
%                         
%                         numPillars = numPillars +1;
%                         SPM(j,6,i) = numPillars; %assign current object new pillar ID
%                         current = SPM(j,9,i); 
%                         [previousIdx,~] = find(SPM(:,6,(i+n)-SPM(SPM(j,4,i),7,i+n))==SPM(SPM(j,4,i),6,i+n));
%                         previous = SPM(previousIdx(1,1),9,(i+n)-SPM(SPM(j,4,i),7,i+n));
%                         if current > previous                            
%                             SPM(SPM(j,4,i),6,i+n) = SPM(j,6,i);
%                             SPM(SPM(j,4,i),7,i+n) = n;
%                             SPM(SPM(j,4,i),9,i+n) = SPM(j,9,i)+1;
%                             SPM(SPM(j,4,i),10,i+n) = SPM(j,10,i);%store starting xy
%                             SPM(SPM(j,4,i),11,i+n) = SPM(j,11,i);
%                         end
                    end
                end
            end
        end
    end
end
end