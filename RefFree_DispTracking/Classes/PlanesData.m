classdef PlanesData
    properties
        nbors
        preplanes % to speed up plane growth
        raw
        refined
        final
        groups
        loc
        locMean
        locComma
        locFilt
        locFiltList
        locTxt
    end
    methods
        function obj = PlanesData(raw3D)
            try
                %this should separate data into 'Preplanes' to narrow feature
                %lists in growPlanes step
                bins = 0:.8:max(raw3D.Z)+1;
                hcs = histcounts(raw3D.Z,bins);
                centers = find(hcs>prctile(hcs,80));
                count = 0;
                
                %Find Dense regions through Z
                for i = 1:size(centers,2);
                    j = i-count;
                    if j<size(centers,2)
                        if centers(1,j+1) - centers(1,j) == 1
                            centers(1,j) = mean(centers(1,j:j+1));
                            centers(:,j+1) = [];
                            j;
                            count = count+1;
                        end
                    end
                end
                
                for i = 1:size(centers,2)-1
                    ceSpaces(1,i) = mean(centers(1,i:i+1));
                end
                
                %Establish Pre-Plane Limits
                for i = 1:size(centers,2)
                    if i == 1
                        ppRegions(1,1) = 0;
                        ppRegions(1,2) = ceSpaces(1,1) + .5;
                    elseif i == size(centers,2)
                        ppRegions(i,1) = ceSpaces(1,i-1) - .5;
                        ppRegions(i,2) = (max(raw3D.Z)+1)/.8;
                    else
                        ppRegions(i,1) = ceSpaces(1,i-1) - .5;
                        ppRegions(i,2) = ceSpaces(1,i) + .5;
                    end
                end
                ppRegions = ppRegions*.8; %convert index back to microns
                
                %Populate Pre-Planes
                index = 1:size(raw3D.Z)';
                for i = 1:size(ppRegions,1)
                    temp = index(raw3D.Z>ppRegions(i,1) & raw3D.Z<ppRegions(i,2))';
                    obj.preplanes(1:size(temp),i) = temp;
                end
            catch
                
                obj.preplanes(:,1) = 1:max(size(raw3D.Z))';
            end
        end
        %%
        
        function obj = nborsPlanes(obj,raw3D,radXY,radZ)
            
            clear obj.nbors
            disp('Finding Neighbors')
            radXY = radXY;
            for i = 1:size(raw3D.X,1)
                %iteratitvely finds all markers in close proximity to
                %marker(i) using input search window
                topX = raw3D.X(i)+ radXY;
                botX = raw3D.X(i)- radXY;
                topY = raw3D.Y(i)+ radXY;
                botY = raw3D.Y(i)- radXY;
                topZ = raw3D.Z(i)+ radZ;
                botZ = raw3D.Z(i)- radZ;
                temp = find(raw3D.X(:)<topX & raw3D.X(:)>botX & raw3D.Y(:)<topY & raw3D.Y(:)>botY& raw3D.Z(:)<topZ & raw3D.Z(:)>botZ);
                obj.nbors(i,1:size(temp)) = temp;
                
            end
            
        end
        
        
        %%
        function obj = nborsPlanesF(obj,raw3D,radXY,radZ)
            %%
            clear obj.nbors
            disp('Finding Neighbors')
            radXY = radXY*2;
            for j = 1:size(obj.preplanes,2)
                pptemp = obj.preplanes(1:nnz(obj.preplanes(:,j)),j);
                dTemp = raw3D.r(pptemp,1:3);
                for i = 1:size(pptemp,1)
                    
                    clear temp
                    %iteratitvely finds all markers in close proximity to
                    %marker(i) using input search window
                    topX = dTemp(i,1)+ radXY;
                    botX = dTemp(i,1)- radXY;
                    topY = dTemp(i,2)+ radXY;
                    botY = dTemp(i,2)- radXY;
                    topZ = dTemp(i,3)+ radZ;
                    botZ = dTemp(i,3)- radZ;
                    
                    temp = pptemp(dTemp(:,1)<topX & dTemp(:,1)>botX & dTemp(:,2)<topY & dTemp(:,2)>botY& dTemp(:,3)<topZ & dTemp(:,3)>botZ);
                    
                    if size(obj.nbors,1) < pptemp(i)
                        obj.nbors(pptemp(i),1:size(temp)) = temp;
                    elseif nnz(obj.nbors(pptemp(i),:)) == 0
                        obj.nbors(pptemp(i),1:size(temp)) = temp;
                    else
                        tempidx = nnz(obj.nbors(pptemp(i),:))+1;
                        tempidx2 = (tempidx+size(temp,1))-1;
                        obj.nbors(pptemp(i),tempidx:tempidx2) = temp;
                    end
                    
                end
            end
            
        end
        
        
        
        function obj = growPlanes(obj,raw3D)
            clear planesTemp
            disp('Growing Planes')
            
            working = 1;
            searched = 1:raw3D.l;
            ss =size(searched,1);
            %start at first row in preplane
            planes = obj.nbors(1,1:nnz(obj.nbors(1,:)))';
            progressbar('Growing Planes')
            j=1; %designates starting at plane 1
            
            while working == 1
                progressbar(size(find(planes(:)),1)/ss)
                newlist = intersect(searched,planes);
                if size(newlist,2)>0
                    for i = 1:size(newlist,1)
                        if i == 1
                            clear new
                            searched((newlist(i,1)==searched)) = [];
                            new(:,1) = obj.nbors(newlist(i,1),1:nnz(obj.nbors(newlist(i,1),:)));
                            
                        else
                            clear newtemp
                            searched((newlist(i,1)==searched)) = [];
                            newtemp(:,1) = obj.nbors(newlist(i,1),1:nnz(obj.nbors(newlist(i,1),:)));
                            new = cat(1,new,newtemp);
                        end
                    end
                end
                sBefore = size(planes,1);
                planes = unique(cat(1,planes,new));
                sAfter = size(planes,1);
                if sBefore == sAfter
                    %if j == 1
                    obj.raw(1:size(planes,1),j) = planes(:,1);
                    j=j+1;
                    %                     elseif size(intersect(obj.raw(:,j-1),planes),2)>0
                    %                             intersect(obj.raw(:,j-1),planes)
                    %                             planes = unique(cat(1,planes(:,1),obj.raw(:,j-1)));
                    %                             obj.raw(1:size(planes,1),j-1) = planes(:,1);
                    %                     else
                    %                     obj.raw(1:size(planes,1),j) = planes(:,1);
                    %                     j=j+1;
                    %                     end
                    
                    clear planes
                    %check to see if any unmatched objects exist
                    
                    finCheck = find(searched,1,'first');
                    if size(finCheck,2) ==0
                        working = 0;
                    else
                        k = searched(find(searched,1,'first'));
                        planes(:,1) = obj.nbors(k,1:nnz(obj.nbors(k,:)));
                        
                    end
                end
            end
        end
        
        
        function [obj,r] = cleanPlanes(obj,raw3D)
            
            %The goal of this function is to remove objects that are
            %members of planes with very few members. This should clean
            %outputs be removing "noise" or unreliable detections.
            r = raw3D.r;
            count =1;
           obj.final = obj.raw;
            
            [planesLoc2,planesGroups,planeSizes] = updateSizes(obj,r);            
            
            %-----------------
            %Remove planes with less than 50 members
            for i = 1:size(planesGroups,1)
                if i <= size(planesGroups,1)
                    planeIdx = planesGroups(i,1:nnz(planesGroups(i,:))); % assign group member planes a number
                    %if current group has less than 50 members, delete all members.
                    for j = 1:nnz(planeIdx)
                        if planeSizes(planeIdx(j)) < 50
                            for k = 1:nnz(obj.raw(:,planeIdx(j)))
                                r(obj.final(k,planeIdx(j)),:) =[]; %clear row in r                                
                                obj.final((obj.final>obj.final(k,planeIdx(j)))) = obj.final((obj.final>obj.final(k,planeIdx(j))))-1; % repeat for obj.final
                            end
                            remove(count) = planeIdx(j);
                            count = count +1;
                        end
                    end
                end                
            end
            try
            obj.final(:,remove) = [];
            catch
            end
            [planesLoc2,planesGroups,planeSizes] = updateSizes(obj,r)
            count = 1;
            clear remove
            
            %-----------------
            %Remove top or bottom plane if there are too few members
            for i = 1:size(planesGroups,1)
                if i <= size(planesGroups,1)
                    planeIdx = planesGroups(i,1:nnz(planesGroups(i,:))); % assign group member planes a number
                    %if current grouped planes are (either the top or bottom
                    %plane, and have less than half the members of the largest
                    %planes), delete all members.
                    if  (nnz(obj.final(:,planeIdx))<(max(planeSizes))/2) && (planesLoc2(i) == max(planesLoc2) || planesLoc2(i) == min(planesLoc2))
                        for j = 1:nnz(planeIdx)
                            for k = 1:nnz(obj.final(:,planeIdx(j)))
                                r(obj.final(k,planeIdx(j)),:) =[];                                
                                obj.final((obj.final>obj.final(k,planeIdx(j)))) = obj.final((obj.final>obj.final(k,planeIdx(j))))-1;
                                remove(count) = planeIdx(j);
                                count = count +1;
                            end
                        end
                    end
                end                
            end
            try
            obj.final(:,remove) = [];
            catch
            end
            [planesLoc2,planesGroups,planeSizes] = updateSizes(obj,r)
            %----------------
            
            obj.final = unique(obj.final','rows')';
            
            function [planesLoc2,planesGroups,planeSizes] = updateSizes(obj,r)
                
                for m = 1:size(obj.final,2)
                    planeSizes(m) = nnz(obj.final(:,m));
                end
                
                %Determine average plane location
                %size(obj.raw,2)
                for m = 1:size(obj.final,2)
                    clear temp
                    temp = obj.final(1:nnz(obj.final(:,m)),m);
                    planesLoc(m) = mean(mean(r(temp,3)));
                end
                
                %Group planes which are close in Z. This should alleviate
                %discrepencies occuring due to cells spanning two or more
                %sections of patterned regions.
                clear planesGroups
                for m = 1:size(planesLoc,2)
                    clear differences
                    differences = planesLoc - planesLoc(1,m);
                    planesGroups(m,1:size(find(abs(differences)<2),2)) = find(abs(differences)<2)';
                end
                planesGroups = unique(planesGroups,'rows');
                
                for m = 1:size(planesGroups)
                    planesGroupsSize(m,1) = nnz(planesGroups(m,:));
                end
                %find planes that are members of more than 1 group and remove
                %them from the smaller groups
                Ui = unique(planesGroups(:));
                U = Ui(1<histc(planesGroups(:),unique(planesGroups(:))));
                for m = 1:size(U,1)
                    clear U2 bigGroup
                    for n = 1:size(planesGroups,1)
                        U2(n) = ismember(U(m),planesGroups(n,:));
                    end
                    if nnz(U2)>0
                        bigGroup = max(planesGroupsSize(U2));
                        
                        for n = 1:size(planesGroups,1)
                            if planesGroupsSize(n,1) ~= bigGroup && U2(n) == 1
                                clear tempGroup
                                tempGroup = planesGroups(n,1:planesGroupsSize(n,1));
                                tempGroup(tempGroup==U(m)) = [];
                                tempGroup(1,planesGroupsSize(n,1)) = 0;
                                planesGroups(n,1:planesGroupsSize(n,1)) = tempGroup(1,1:planesGroupsSize(n,1));
                            end
                        end
                    end
                end
                
                %Repeat average Z plane location with merged planes
                for m = 1:size(planesGroups,1)
                    planesLoc2(m) = mean(planesLoc(planesGroups(m,1:nnz(planesGroups(m,:)))));
                end
            end
        end
    end
end