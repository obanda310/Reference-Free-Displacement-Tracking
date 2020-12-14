function [lub] = refreshData(lub,maxAi,maxAc,numPillars,maxMasks)
borders = double(zeros(size(maxMasks,1),size(maxMasks,2)));
borders(1:10,:) = 1;
borders(:,1:10) = 1;
borders(end-9:end,:) = 1;
borders(:,end-9:end) = 1;
pillarCO =(max(lub(:,6)/2));

for i = 1:numPillars
    pillarStart = find(lub(:,7)==i,1,'first');
    pillarEnd = find(lub(:,7)==i,1,'last');
    pillarSize = ((pillarEnd+1)-pillarStart);% # of objects in Pillar
    pillarSize2 = lub(pillarEnd,6)-lub(pillarStart,6);% Height of Pillar
    lub(pillarStart:pillarEnd,27) = pillarSize;
    lub(pillarStart:pillarEnd,21) = pillarSize2; % Height of Pillar
    lub(pillarStart:pillarEnd,9) = 1:1:pillarSize; %position within pillar (not necessarily same as frame number)
    lub(pillarStart:pillarEnd,10) = lub(pillarStart,1); %pillar x start
    lub(pillarStart:pillarEnd,11) = lub(pillarStart,2); %pillar y start
    lub(pillarStart:pillarEnd,24) = lub(pillarStart,10);
    lub(pillarStart:pillarEnd,25) = lub(pillarStart,11);
    for j = 2:pillarSize
        lub(pillarStart+j-1,24) = mean(lub(pillarStart:pillarStart+j-1,1));
        lub(pillarStart+j-1,25) = mean(lub(pillarStart:pillarStart+j-1,2));
    end
end
pillarSizeLimit = round(max(lub(:,21))/2); %minimum pillar size deemed acceptable
lub(:,12) = sqrt((lub(:,1)-lub(:,10)).^2 +(lub(:,2)-lub(:,11)).^2); %Distance from start

% Identifying problem Pillars
for i = 3:size(lub,1)
    lub(i,13) = sqrt((lub(i,1)-lub(i-1,1))^2 +(lub(i,2)-lub(i-1,2))^2); %Distance from previous
    
    %This code should determine if trajectory from the start to the current
    %is reasonable.
    if (lub(i,12) > 5) && (lub(i,9) > 1) && ((lub(i,13) > 2) || lub(i,6)>pillarCO)
        lub(i,14) = atan2d(norm(cross([lub(i,1:2),0]-[lub(i,10:11),0],[lub(i,24:25),0]-[lub(i,10:11),0])),dot([lub(i,1:2),0]-[lub(i,10:11),0],[lub(i,24:25),0]-[lub(i,10:11),0])); %Total Angle to [lub(i,1:2),0]
        if lub(i,14)> maxAi
            lub(i,15) = 1; %If total angle is too large, +1 to score
        end
    else %Do not calculate angles for brownian movement
        lub(i,14) = 0;
        lub(i,15) = 0;
    end
    
    %This code should determine if the most recent trajectory is reasonable
    if ((lub(i,13) >3) || (lub(i,6)>pillarCO && lub(i,13) >3)) && (lub(i,6) ~=1) && (lub(i,9) > 2)
        lub(i,16) = dot([lub(i,24:25),0]-[lub(i,10:11),0],[lub(i,1:2),0]-[lub(i-1,1:2),0]); %Relationship of current traj to traj history (if negative, object is moving in wrong direction!)
        lub(i,18) = atan2d(norm(cross([lub(i,1:2),0]-[lub(i-2,1:2),0],[lub(i-1,1:2),0]-[lub(i-2,1:2),0])),dot([lub(i,1:2),0]-[lub(i-2,1:2),0],[lub(i-1,1:2),0]-[lub(i-2,10:11),0])); %pillar end angle
        if lub(i,16)<0
            lub(i,17) = 1; %If trajectory is negative (approaching opposite of previous direction), flag it
        elseif lub(i,16)>0 && lub(i,18)>maxAc
            lub(i,19) = 1; %If ending angle is too large, flag it
        end
    else
        lub(i,16) = 0;
        lub(i,17) = 0;
        lub(i,18) = 0;
        lub(i,19) = 0;
    end
    
    
    %If a pillar is too small, flag it
    if lub(i,21) < pillarSizeLimit || lub(i,27)<3
        lub(i,22) = 1;
    else
        lub(i,22) = 0;
    end
    
    %If a pillar is near the edge, ignore it
    if borders(round(lub(i,2)),round(lub(i,1)))==1
        lub(i,26) = 0;
    else
        lub(i,26) = 1;
    end
    
    %Calculate a problem score
    lub(i,20) = (lub(i,22))*lub(i,26); %total Score  lub(i,15)+ lub(i,17)+lub(i,19)+
end



maxPillarCO = max(lub(:,21))-2;
for i = 1:numPillars
    clear tempPillar
    tempPillar = lub(lub(:,7)==i,:);
    try
        if sum(tempPillar(:,20)) <= 0 && max(tempPillar(:,12))<5 && max(tempPillar(:,21))>maxPillarCO
            lub(lub(:,7)==i,30) = 1;
        end
    catch
    end
    
end

end