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
            %start at first row in preplane
            planes = obj.nbors(1,1:nnz(obj.nbors(1,:)))';
            
            j=1; %designates starting at plane 1
            while working == 1
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
            clear planesFinal
            r = raw3D.r;
            j =1;
            obj.final = 0;
            for i = 1:size(obj.raw,2)
                if nnz(obj.raw(:,i))>50
                    obj.final(1:nnz(obj.raw(:,i)),j) = obj.raw(1:nnz(obj.raw(:,i)),i);
                    j=j+1;
                else
                    for k = 1:nnz(obj.raw(:,i))
                        r(obj.raw(k,i),:) =[];
                        obj.raw((obj.raw>obj.raw(k,i))) = obj.raw((obj.raw>obj.raw(k,i)))-1;
                        obj.final((obj.final>obj.raw(k,i))) = obj.final((obj.final>obj.raw(k,i)))-1;
                    end
                end
            end
        end
    end
end