function [subpixMaxima,noPillars] = LinkMaxima(subpixMaxima,maxLD,maxJD)
% Find number of objects per frame in subpixMaxima
for i = 1:size(subpixMaxima,3)
    [lastNZElement,~] = find(subpixMaxima(:,1,i),1,'last');
    twoDimMaxIndSub(i,1) = lastNZElement;
end

%closest neighbor in frame above
noPillars = 0;
for i = 1:size(subpixMaxima,3)
    for j = 1:(twoDimMaxIndSub(i,1))
        
        %if on first frame, assign the current object a new pillar index
        if i == 1
            noPillars = noPillars +1;
            subpixMaxima(j,6,i) = noPillars;
            subpixMaxima(j,9,i) = 1;
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
                    
                    %if a pillar exists for current object but not it's
                    %nearest neighbor in following frame, assign it's
                    %nearest neighbor the pillar
                    if subpixMaxima(j,6,i) > 0 && subpixMaxima(subpixMaxima(j,4,i),6,i+n) == 0
                        subpixMaxima(subpixMaxima(j,4,i),6,i+n) = subpixMaxima(j,6,i);
                        subpixMaxima(subpixMaxima(j,4,i),7,i+n) = n;
                        subpixMaxima(subpixMaxima(j,4,i),9,i+n) = subpixMaxima(j,9,i)+1;
                        
                        %if a pillar exists for current object and also it's
                        %nearest neighbor in following frame, determine
                        %whether current pillar is a better match. If it
                        %is, replace the nearest neighbor pillar
                    elseif subpixMaxima(j,6,i) > 0 && subpixMaxima(subpixMaxima(j,4,i),6,i+n) > 0
                        current = subpixMaxima(j,9,i);
                        [previousIdx,~] = find(subpixMaxima(:,6,(i+n)-subpixMaxima(subpixMaxima(j,4,i),7,i+n))==subpixMaxima(subpixMaxima(j,4,i),6,i+n));
                        previous = subpixMaxima(previousIdx(1,1),9,(i+n)-subpixMaxima(subpixMaxima(j,4,i),7,i+n));
                        if current > previous
                            subpixMaxima(subpixMaxima(j,4,i),6,i+n) = subpixMaxima(j,6,i);
                            subpixMaxima(subpixMaxima(j,4,i),7,i+n) = n;
                            subpixMaxima(subpixMaxima(j,4,i),9,i+n) = subpixMaxima(j,9,i)+1;
                        end
                        
                        %if no pillar exists for current object or it's
                        %nearest neighbor in following frame, assign it
                        %a new pillar and assign it's nearest neighbor the
                        %pillar as well
                    elseif subpixMaxima(j,6,i) == 0 && subpixMaxima(subpixMaxima(j,4,i),6,i+n) == 0
                        noPillars = noPillars +1;
                        subpixMaxima(j,6,i) = noPillars;
                        subpixMaxima(subpixMaxima(j,4,i),6,i+n) = subpixMaxima(j,6,i);
                        subpixMaxima(subpixMaxima(j,4,i),7,i+n) = n;
                        subpixMaxima(j,9,i) = 1;
                        subpixMaxima(subpixMaxima(j,4,i),9,i+n) = subpixMaxima(j,9,i)+1;
                        
                        %if no pillar exists for current object, but one
                        %does exist for it's nearest neighbor in following
                        %frame, create new pillar, and determine if current
                        %pillar is a better match
                    elseif subpixMaxima(j,6,i) == 0 && subpixMaxima(subpixMaxima(j,4,i),6,i+n) > 0
                        noPillars = noPillars +1;
                        subpixMaxima(j,6,i) = noPillars;
                        current = subpixMaxima(j,9,i);
                        [previousIdx,~] = find(subpixMaxima(:,6,(i+n)-subpixMaxima(subpixMaxima(j,4,i),7,i+n))==subpixMaxima(subpixMaxima(j,4,i),6,i+n));
                        previous = subpixMaxima(previousIdx(1,1),9,(i+n)-subpixMaxima(subpixMaxima(j,4,i),7,i+n));
                        if current > previous
                            subpixMaxima(subpixMaxima(j,4,i),6,i+n) = subpixMaxima(j,6,i);
                            subpixMaxima(subpixMaxima(j,4,i),7,i+n) = n;
                            subpixMaxima(subpixMaxima(j,4,i),9,i+n) = subpixMaxima(j,9,i)+1;
                        end
                    end
                end                
            end
        end               
    end
end
end