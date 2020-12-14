function [SPM,numPillars] = LinkMaxima(SPM,maxLD,maxJD)
count = 0;


maxAi = 20; %Maximum angle deviation measured from initial location and immediately previous location
maxAc = 30; %Maximum angle deviation measured from two immediately previous locations
maxTD = maxLD*2;
maxS = 3; %Maximum deviance score allowable for a match


%SPM is subpixMaxima from tracking code
%LastNZ = 2-D Maximum Indices of SPM (number of objects on a page of SPM)

%SPM Structure
%(x,:,:) describes object ID
%(:,x,:) describes object parameter (the indices below)
%(:,:,x) describes object frame

%SPM Indices
%1 X
%2 Y
%3 Z (frame)
%4 Index of Following Match
%5 Number of Frames Up to Following Match
%6 Pillar ID
%7 Number of Frames Down to Previous Match
%8 Index of Previous Match
%9 Current objects place in pillar (to measure the number of members in a
%pillar)
%10 Starting X in pillar
%11 Starting Y in pillar
%12 Angle in degrees (ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));)
%13 Competitor's angle in degrees
%14 XY distance to nearest neighbor above
%15 Match Score Foreward
%16 Match Score Backward

% Find number of objects per frame in SPM
for i = 1:size(SPM,3)
    [lastNZElement,~] = find(SPM(:,1,i),1,'last');
    LastNZ(i,1) = lastNZElement;
end

%%%%%%%%%%%%%%%%%%%%%PILLAR LINKING START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up Outer Looping Sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numPillars = 0;
for i = 1:size(SPM,3)
    count1 = 0;
    count2 = 0;
    count3 = 0;
    count4 = 0;
    count5 = 0;
    A = strcat('Linking Frame:  ',num2str(i),'****************');
    disp(A)
    for j = 1:(LastNZ(i,1))
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Setting up First Pillars on First Frame
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %if on first frame, assign the current object a new pillar index
        %and continue
        if i == 1
            numPillars = numPillars +1; %increase pillar count
            SPM(j,6,i) = numPillars; %store pillar ID
            SPM(j,9,i) = 1; %record current pillar length (has to be 1 on first frame)
            SPM(j,10,i) = SPM(j,1,i); %store first frame of pillar XY --
            SPM(j,11,i) = SPM(j,2,i); %--redundant, but storing for consistency later
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Setting up Inner Looping Sequence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ignoreN = 0; %end the loop early if this value changes
        for n = 1:maxJD %check all frames within jump distance range
            
            if ignoreN ==0 && i+n <= size(SPM,3)
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Establish the preliminary list of best match in following
                %frame
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                clear tempDistances nearUpNeighbors nearUpNeighbors2 nearUpNeighborFinal
                %calulate distances from maxima in frame above to current object
                tempDistances = sqrt((SPM(1:LastNZ(i,1),1,i+n)-SPM(j,1,i)).^2 +(SPM(1:LastNZ(i,1),2,i+n)-SPM(j,2,i)).^2);
                
                if SPM(j,10,i)>0
                    currentDistance = sqrt((SPM(j,1,i)-SPM(j,10,i))^2 +(SPM(j,2,i)-SPM(j,11,i))^2);
                    %                     if currentDistance > maxLD
                    %                         disp(num2str(currentDistance))
                    %                     end
                else
                    currentDistance = 1;
                end
                
                %name the object in the next frame with the least distance from
                %current object its "nearest neighbor"
                if currentDistance > 3
                    %find all in larger range
                    [nearUpNeighbors2,~] = find(tempDistances<maxTD);
                    
                    %create triangular prediction of prospective match
                    %location
                    dScale = maxTD/currentDistance;
                    mid(1,1) = SPM(j,1,i)+dScale*(SPM(j,1,i)-SPM(j,10,i));
                    mid(1,2) = SPM(j,2,i)+dScale*(SPM(j,2,i)-SPM(j,11,i));
                    perpD = (sqrt((mid(1,1)-SPM(j,10,i))^2 +(mid(1,2)-SPM(j,11,i))^2))*(tan(((maxAi/2))*(pi/180)));                   
                    t = sqrt((mid(1,1)-SPM(j,10,i))^2+(mid(1,2)-SPM(j,11,i))^2); %or currentDistance + maxTD
                    Tri(3,1) = mid(1,1)-((mid(1,2)-SPM(j,11,i))/t)*perpD;
                    Tri(3,2) = mid(1,2)+((mid(1,1)-SPM(j,10,i))/t)*perpD;
                    Tri(1,1) = mid(1,1)+((mid(1,2)-SPM(j,11,i))/t)*perpD;
                    Tri(1,2) = mid(1,2)-((mid(1,1)-SPM(j,10,i))/t)*perpD;
                    Tri(2,1) = SPM(j,10,i);
                    Tri(2,2) = SPM(j,11,i);
                    plot(Tri(:,1),Tri(:,2))
                    hold on
                    %filter nearUpNeighbors to include only objects within
                    %the prediction triangle
                    
                    in = inpolygon(SPM(nearUpNeighbors2,1,i+n),SPM(nearUpNeighbors2,2,i+n),Tri(:,1),Tri(:,2));
                    if in == 0
                    [nearUpNeighbors,~] = find(tempDistances<maxLD);
                    else
                    [nearUpNeighbors3,~] = find(tempDistances<maxLD);
                    nearUpNeighbors = cat(1,nearUpNeighbors2(in),nearUpNeighbors3);
                    end
                    
                else
                    [nearUpNeighbors,~] = find(tempDistances<maxLD);
                end
                
                nearUpNeighbors(:,2) = tempDistances(nearUpNeighbors);
                
                AiW = (i+(n-1))/size(SPM,3); % A weight factor for total angles
                %based on current frame (upper frames should have most of
                %the deformation, so mute the angle influences until later
                %frames)
                
                
                clear AcW
                AcW = ((currentDistance/maxTD)); % A weight factor for current angles
                %based on total distance traveled (persistence is measured
                %best when there is signficant distance traveled away from
                %original location)
                %                 if AcW>1
                %                     disp(num2str(AcW))
                %                     disp(num2str(SPM(j,10,i)))
                %                 end
                
                
                
                
                if nearUpNeighbors > 0
                    skip =0;
                    if SPM(j,9,i) > 2 %pillar length must be greater than 2 to have enough points to create an angle
                        noScore = 0;
                        for m = 1:size(nearUpNeighbors)
                            Ini = [SPM(j,10:11,i) 0]; %Corner - First in Pillar
                            Prev = [SPM(SPM(j,8,i),1:2,i-SPM(j,7,i)) 0]; %Corner - Previous Object in Pillar
                            Curr = [SPM(j,1:2,i) 0]; %Current
                            Pros = [SPM(nearUpNeighbors(m,1),1:2,i+n) 0]; %Prospective
                            nearUpNeighbors(m,3) = sqrt((SPM(nearUpNeighbors(m,1),1,i+n)-SPM(j,10,i)).^2 +(SPM(nearUpNeighbors(m,1),2,i+n)-SPM(j,11,i)).^2)/maxTD;
                            nearUpNeighbors(m,4) = atan2d(norm(cross(Pros-Ini,Prev-Ini)),dot(Pros-Ini,Prev-Ini)); %Total Angle to Pros
                            nearUpNeighbors(m,5) = atan2d(norm(cross(Pros-Prev,Curr-Prev)),dot(Pros-Prev,Curr-Prev)); %Pillar-end Angle
                            nearUpNeighbors(m,6) = atan2d(norm(cross(Curr-Ini,Prev-Ini)),dot(Curr-Ini,Prev-Ini)); %Total Angle to Current
                            
                            direction = dot(Curr-Ini,Pros-Curr);
                            if direction < 0 && currentDistance > 2 && nearUpNeighbors(m,2) > 1
                                dCheck = 2;
                            else
                                dCheck = 0;
                            end
                            
                            %Check for persistence in trajectory ~ if there
                            %is a history i.e. the object has moved from
                            %its original location. Small values should
                            %imply that the object persists in one
                            %direction.
                            if currentDistance > 2 && nearUpNeighbors(m,4) < maxAi
                                trajCheck = 0;
                                %trajCheck = (round((ceil(nearUpNeighbors(m,4))/ceil(nearUpNeighbors(m,6)))));
                                %disp(num2str(trajCheck))
                            elseif currentDistance > 2 && nearUpNeighbors(m,4) > maxAi
                                trajCheck = 1;
                            else
                                trajCheck = 0;
                            end
                            
                            
                            if nearUpNeighbors(m,2) > maxLD/2
                                distCheck = 0.5;
                            else
                                distCheck = 0;
                            end
                            
                            
                            %Handles abrupt changes in trajectory of
                            %objects with a history
                            if nearUpNeighbors(m,2) > .3 && currentDistance >2 && nearUpNeighbors(m,5) > maxAc
                                noiseCheck = 2;
                            %Handles abrupt changes in location and
                            %trajectory in a previously stable pillar
                            elseif nearUpNeighbors(m,2) > 2 && currentDistance < 2 && nearUpNeighbors(m,4) > maxAi
                                noiseCheck = 2;
                            else
                            %If no abrupt change, keep it at zero!    
                                noiseCheck = 0;
                            end
                            
                            %Benefits matches under a threshold distance
                            if nearUpNeighbors(m,2)<.5
                                closeCheck = -1;
                            else
                                closeCheck = 0;
                            end
                            
                            %Penalizes Matches exceeding a certain total
                            %distance
                            if nearUpNeighbors(m,3)>(maxLD*2)
                                limitCheck = 3;
                            else
                                limitCheck = 0;
                            end
                            
                            %Create quality scores here:
                            
                            nearUpNeighbors(m,7) = dCheck + noiseCheck + trajCheck + distCheck +closeCheck +limitCheck;
               
                        end
                        
                        nearUpNeighbors = sortrows(nearUpNeighbors,[7 2]);
                        nearUpNeighborFinal = nearUpNeighbors((nearUpNeighbors(:,7)<=maxS),1);
                        %                         if size(nearUpNeighborFinal,1) >1
                        %                             disp(nearUpNeighborFinal)
                        %                             disp(num2str(i))
                        %                         end
                        if min(nearUpNeighbors(:,7))>maxS
                            skip =1;
                            %disp(num2str(min(nearUpNeighbors(:,6))))
                            %disp(num2str(size(nearUpNeighbors(:,6),1)))
                            %disp('Skipped')
                        end
                    else
                        %the simple case (not a lot of motion, no pillar history, so no angles)
                        noScore = 1;
                        nearUpNeighborFinal = nearUpNeighbors((nearUpNeighbors(:,2)==min(nearUpNeighbors(:,2))),1);
                        count5 = count5 + 1;
                    end
                else
                    skip = 1;
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Pillar Assignment and Checks Against Previous Pillar Match Assignments
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if skip == 0
                    %if the nearest neighbor falls within the max linking distance,
                    if size(nearUpNeighborFinal,1)>0
                        matched = 0;
                        for o = 1:size(nearUpNeighborFinal)
                            if matched == 0
                                %store the nearest up neighbor as as the fourth index in
                                %SPM
                                SPM(j,4,i) = nearUpNeighborFinal(o,1);
                                SPM(j,5,i) = n;
                                if noScore == 0
                                    SPM(j,15,i) = nearUpNeighbors(o,7);
                                else
                                    %New pillars created late in the game
                                    %deserve a bad score!
                                    SPM(j,15,i) = 1000;
                                end
                                %store the distance to the nearest up
                                %neighbor as index 14
                                SPM(j,14,i) = tempDistances(nearUpNeighborFinal(o,1),1);
                                
                                
                                
                                
                                %if a pillar exists for current object and also it's
                                %nearest neighbor in following frame, determine
                                %whether current pillar is a better match. If it
                                %is, replace the nearest neighbor pillar
                                if SPM(j,6,i) > 0 && SPM(SPM(j,4,i),6,i+n) > 0
                                    count1 = count1 +1;
                                    %count = count+1;
                                    %disp(num2str(count))
                                    %Compare Pillar Size First! (don't
                                    %compete with small pillars!)
                                    if SPM(j,9,i) > (2 * SPM(SPM(j,4,i),9,i+n))
                                        SPM(SPM(j,4,i),6,i+n) = SPM(j,6,i); %assign current pillar ID
                                        SPM(SPM(j,4,i),7,i+n) = n; %store n to find nearest neighbor frame later if needed
                                        SPM(SPM(j,4,i),8,i+n) = j; %store i to find nearest neighbor frame later if needed
                                        SPM(SPM(j,4,i),9,i+n) = SPM(j,9,i)+1; %store pillar size as current size +1
                                        SPM(SPM(j,4,i),10,i+n) = SPM(j,10,i);%store starting xy
                                        SPM(SPM(j,4,i),11,i+n) = SPM(j,11,i);
                                        if noScore == 0
                                            SPM(SPM(j,4,i),16,i) = min(nearUpNeighbors(:,7));
                                        else
                                            SPM(SPM(j,4,i),16,i) = 10;
                                        end
                                        matched = 1;
                                        ignoreN = 1; %a match has been found the loop should proceed to next object after pillar assignment                                    
                                    
                                    % Next, compare current match score with previous match score
                                    
                                    elseif SPM(j,15,i) < SPM(SPM(j,4,i),16,i+n)
                                        
                                        SPM(SPM(j,4,i),6,i+n) = SPM(j,6,i); %assign current pillar ID
                                        SPM(SPM(j,4,i),7,i+n) = n; %store n to find nearest neighbor frame later if needed
                                        SPM(SPM(j,4,i),8,i+n) = j; %store i to find nearest neighbor frame later if needed
                                        SPM(SPM(j,4,i),9,i+n) = SPM(j,9,i)+1; %store pillar size as current size +1
                                        SPM(SPM(j,4,i),10,i+n) = SPM(j,10,i);%store starting xy
                                        SPM(SPM(j,4,i),11,i+n) = SPM(j,11,i);
                                        if noScore == 0
                                            SPM(SPM(j,4,i),16,i) = min(nearUpNeighbors(:,7));
                                        else
                                            SPM(SPM(j,4,i),16,i) = 10;
                                        end
                                        matched = 1;
                                        ignoreN = 1; %a match has been found the loop should proceed to next object after pillar assignment
                                    else
                                        count4 = count4 +1;
                                        matched = 0;
                                        ignoreN = 0; %a match has been found the loop should proceed to next object after pillar assignment
                                    end
                                    
                                    %Finally, select another match, if
                                    %available, for the pillar whose match
                                    %was replaced, and remove the replaced
                                    %match from its list of candidates
                                    
                                    
                                    
                                    
                                    %if a pillar exists for current object but not it's
                                    %nearest neighbor in following frame, assign it's
                                    %nearest neighbor the pillar
                                elseif SPM(j,6,i) > 0 && SPM(SPM(j,4,i),6,i+n) == 0
                                    count2 = count2 +1;
                                    SPM(SPM(j,4,i),6,i+n) = SPM(j,6,i); %assign current pillar ID
                                    SPM(SPM(j,4,i),7,i+n) = n; %store n to find nearest neighbor frame later if needed
                                    SPM(SPM(j,4,i),8,i+n) = j; %store i to find nearest neighbor frame later if needed
                                    SPM(SPM(j,4,i),9,i+n) = SPM(j,9,i)+1; %store pillar size as current size +1
                                    SPM(SPM(j,4,i),10,i+n) = SPM(j,10,i);%store starting xy
                                    SPM(SPM(j,4,i),11,i+n) = SPM(j,11,i);
                                    if noScore == 0
                                        SPM(SPM(j,4,i),16,i) = min(nearUpNeighbors(:,7));
                                    else
                                        SPM(SPM(j,4,i),16,i) = 10;
                                    end
                                    
                                    %if there are enough points to calculate an angle
                                    %ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
                                    if SPM(SPM(j,4,i),9,i+n) > 2
                                        P0 = SPM(SPM(j,4,i),10:11,i+n);
                                        P1 = SPM(j,1:2,i);
                                        P2 = SPM(SPM(j,4,i),1:2,i+n);
                                        SPM(SPM(j,4,i),12,i+n) = (180/pi)*atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
                                    end
                                    
                                    matched = 1;
                                    ignoreN = 1; %a match has been found the loop should proceed to next object after pillar assignment
                                    
                                    %if no pillar exists for current object or it's
                                    %nearest neighbor in following frame, assign it
                                    %a new pillar and assign it's nearest neighbor the
                                    %pillar as well
                                elseif SPM(j,6,i) == 0 && SPM(SPM(j,4,i),6,i+n) == 0
                                    count3 = count3 +1;
                                    numPillars = numPillars +1;
                                    SPM(j,6,i) = numPillars;
                                    SPM(SPM(j,4,i),6,i+n) = numPillars;
                                    SPM(SPM(j,4,i),7,i+n) = n;
                                    SPM(SPM(j,4,i),8,i+n) = j;
                                    SPM(j,9,i) = 1; %current object is pillar size 1
                                    SPM(SPM(j,4,i),9,i+n) = 2; %neighbor object makes pillar size 2
                                    SPM(j,10,i) = SPM(j,1,i);%store starting xy
                                    SPM(j,11,i) = SPM(j,2,i);
                                    SPM(SPM(j,4,i),10,i+n) = SPM(j,1,i);%store starting xy
                                    SPM(SPM(j,4,i),11,i+n) = SPM(j,2,i);
                                    if noScore == 0
                                        SPM(SPM(j,4,i),16,i) = min(nearUpNeighbors(:,7));
                                    else
                                        SPM(SPM(j,4,i),16,i) = 10;
                                    end
                                    matched = 1;
                                    ignoreN = 1; %a match has been found the loop should proceed to next object after pillar assignment
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
    end
    disp(['Low-Info Match: ' num2str(count5)])
    %disp(['Attempted Replacements: ' num2str(count1)])
    %disp(['Failed Replacements: ' num2str(count4)])
    disp(['Replacements: ' num2str(count1-count4)])
    %disp(['Continuations: ' num2str(count2)])
    disp(['Initializations: ' num2str(count3)])
end

end




%% Old version of pillar match replacement
%Compare Angles first
%                         if SPM(SPM(j,4,i),9,i+n) > 2 && SPM(j,9,i) > 2
%                             P0 = SPM(SPM(j,4,i),10:11,i+n);
%                             P1 = SPM(j,1:2,i);
%                             P2 = SPM(SPM(j,4,i),1:2,i+n);
%                             SPM(SPM(j,4,i),13,i+n) = (180/pi)*atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
%                             if SPM(SPM(j,4,i),13,i+n) < SPM(SPM(j,4,i),12,i+n)
%                                 SPM(SPM(j,4,i),6,i+n) = SPM(j,6,i);
%                                 SPM(SPM(j,4,i),7,i+n) = n;
%                                 SPM(SPM(j,4,i),8,i+n) = i;
%                                 SPM(SPM(j,4,i),9,i+n) = SPM(j,9,i)+1;
%                             end
%                         else
%                             %Otherwise use pillar length
%                             current = SPM(j,9,i);
%                             previous = SPM(SPM(SPM(j,4,i),8,i+n),9,(i+n)-(SPM(SPM(j,4,i),7,i+n)));
%                             if current > previous
%                                 SPM(SPM(j,4,i),6,i+n) = SPM(j,6,i);
%                                 SPM(SPM(j,4,i),7,i+n) = n;
%                                 SPM(SPM(j,4,i),8,i+n) = i;
%                                 SPM(SPM(j,4,i),9,i+n) = SPM(j,9,i)+1;
%                             end
%                         end