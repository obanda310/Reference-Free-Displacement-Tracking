%This function file is currently not in use as of 07/12/2017. The code
%referencing this function file has been removed from the main tracking
%script.

%Old code from tracking script:
%Copy-Paste into tracking script to use
% for k = 1:4
% numPillars = max(lub(:,7));
% [lub] = refreshData(lub,maxAi,maxAc,numPillars,maxMasks);
% problems = find(lub(:,20)>0);
% % Fixing Problem Pillars
% sLocs = unique(lub(:,[10 11 7]),'rows','stable');
% exempt = unique(lub((lub(:,21)<(0.8*size(roiImgs,3))),7)); % exempt from being matched to because the pillar is too short
% for i = 1:size(exempt,1)
% sLocs(sLocs(:,3)==exempt(i,1),:) = 0;
% end
% sLocs = unique(sLocs,'rows','stable');
% sLocs(sLocs(:,3)==0,:) = [];
% lub(:,23) = lub(:,7);
% [lub] = fixProblemPillars(lub,maxAi,problems,sLocs,maxD);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%This function attempts to correct "problems" automatically. In its current
%state, however, it seems to create just as many problems as it fixes.
function [lub] = fixProblemPillars(lub,maxAi,problems,sLocs,maxD)
count = 0;
for i = 1:size(problems,1)
    clear distances nbors
    currentPillar = lub(problems(i,1),7);
    firstFrame = lub(find(lub(:,7)==currentPillar,1,'first'),6);
    if firstFrame >= max(lub(:,6))-3
        maxD2 = maxD*.5;
    else
        maxD2 = maxD;
    end
    distances = sqrt((sLocs(:,1)-lub(problems(i,1),1)).^2+(sLocs(:,2)-lub(problems(i,1),2)).^2);
    nbors = sLocs((distances<maxD2),3);
    
    for j = 1:size(nbors,1)
%         if size(find((lub(:,7) == nbors(j,1)) & (lub(:,6) < lub(problems(i,1),6)),1,'last'),1)>0
%             nbors(j,2) = find(lub(:,7) == nbors(j,1) & lub(:,6) < lub(problems(i,1),6),1,'last');
%         else
%             nbors(j,2) = find(lub(:,7) == nbors(j,1) & lub(:,6) <= lub(problems(i,1),6),1,'last');
%         end
        
        
        idx = find((lub(:,7)==nbors(j,1))&(lub(:,6)<firstFrame),1,'last');
        if size(idx,1) == 0
            nbors(j,3) = 180;
        else
            nbors(j,3) = atan2d(norm(cross([lub(problems(i,1),1:2),0]-[lub(idx,10:11),0],[lub(idx,1:2),0]-[lub(idx,10:11),0])),dot([lub(problems(i,1),1:2),0]-[lub(idx,10:11),0],[lub(idx,1:2),0]-[lub(idx,10:11),0])); %Total Angle to [lub(i,1:2),0]
        end
        nbors(j,4) = distances(sLocs(:,3)==nbors(j,1)); %distance from problem to neighbors
        %nbors(j,5) = dot([lub(nbors(j,2),24:25),0]-[lub(nbors(j,2),10:11),0],[lub(problems(i,1),1:2),0]-[lub(nbors(j,2),1:2),0]);%Trajectory check. If negative, then going the wrong way
        %nbors(j,6) = sqrt((lub(problems(i,1),1)-lub(nbors(j,2),1))^2+(lub(problems(i,1),2)-lub(nbors(j,2),2))^2);
      
        nbors(j,7) = lub(find(lub(:,7) == nbors(j,1),1,'first'),30);
    end
    %     [del,~] = find(((nbors(:,5)<0) & (nbors(:,6)>= (maxD/2)))|nbors(j,7)==1);
    
    nbors(nbors(:,7)==1,:) = [];
    
    clear nborsSort
    if size(nbors,1) > 0
        nborsSort = sortrows(nbors,3);
        if lub(problems(i,1),23)>0 && lub(nborsSort(1,1),7) == lub(problems(i,1),23)
            nborsSort(1,:)= [];
        end
        
        
        if size(nborsSort,1) > 0
            if nborsSort(1,3) < maxAi
                
                % if an object already exists at replacement position in pillar,
                % remove the previous object pillar designation
                clear replace
                shift = find((lub(:,7) == lub(problems(i,1),7))&(lub(:,6)>=lub(problems(i,1),6)));
                if size(find(lub(:,7) == nborsSort(1,1) & lub(:,6) >= lub(problems(i,1),6)),1)>0
                    %finalIdx = find(lub(:,7)==nborsSort(1,1),1,'last');
                    replace = find(lub(:,7) == nborsSort(1,1) & lub(:,6) >= lub(problems(i,1),6));
                    
                    if replace>0
                        count = count+1;
                        %finalIdx
                    end
                    %lub(replace:finalIdx,23) = nborsSort(1,1); %store previous desgination to avoid changing back later
                    
                    
                    if size(shift,1) == 1 %this should deal with pillars that have more than one object in the same frame. Typically happens when there is very large deformation (~4microns)
                    else
                    lub(replace,7) = max(lub(:,7))+1; %change pillar designation to new value for all replaced objects
                    lub(shift,7) = lub(find(lub(:,7)==nborsSort(1,1),1,'first'),7);
                    lub(shift,10) = lub(find(lub(:,7)==nborsSort(1,1),1,'first'),10);
                    lub(shift,11) = lub(find(lub(:,7)==nborsSort(1,1),1,'first'),11);
                    sLocs(currentPillar,1:2) = 0;
                    end
                else
                    lub(shift,7) = lub(find(lub(:,7)==nborsSort(1,1),1,'first'),7);
                    lub(shift,10) = lub(find(lub(:,7)==nborsSort(1,1),1,'first'),10);
                    lub(shift,11) = lub(find(lub(:,7)==nborsSort(1,1),1,'first'),11);
                    sLocs(currentPillar,1:2) = 0;
                end
            end
        end
    end
end
lub = sortrows(lub,[7 6]);
count
end