classdef PlanesData
    properties
        nbors
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
        function obj = PlanesData(raw3D,radXY,radZ)
            clear obj.nbors
            for i = 1:size(raw3D.X,1)
                %iteratitvely finds all markers in close proximity to
                %marker(i) using input search window
                topX = raw3D.X(i)+ radXY;
                botX = raw3D.X(i)- radXY;
                topY = raw3D.Y(i)+ radXY;
                botY = raw3D.Y(i)- radXY;
                topZ = raw3D.Z(i)+ radZ;
                botZ = raw3D.Z(i)- radZ;
                obj.nbors(i,1:size(find(raw3D.X(:)<topX & raw3D.X(:)>botX & raw3D.Y(:)<topY & raw3D.Y(:)>botY& raw3D.Z(:)<topZ & raw3D.Z(:)>botZ))) = find(raw3D.X(:)<topX & raw3D.X(:)>botX & raw3D.Y(:)<topY & raw3D.Y(:)>botY& raw3D.Z(:)<topZ & raw3D.Z(:)>botZ);
            end
        end
        
        function obj = growPlanes(obj,raw3D)
            clear planesTemp
            working = 1;
            searched = 1:1:raw3D.l;
            %start at first row in r
            planesTemp(:,1) = obj.nbors(1,1:nnz(obj.nbors(1,:)));
            planes = planesTemp;
            j=1; %designates starting at plane 1
            while working == 1
                for i = 1:size(planes)
                    if ismember(planes(i,1),searched) == 1
                        clear new
                        searched((planes(i,1)==searched)) = [];
                        new(:,1) = obj.nbors(planes(i,1),1:nnz(obj.nbors(planes(i,1),:)));
                        planesTemp = cat(1,planesTemp,new);
                    end
                end
                sBefore = size(planes,1);
                planes = unique(cat(1,planes,planesTemp));
                sAfter = size(planes,1);
                if sBefore == sAfter
                    obj.raw(1:size(planes,1),j) = planes(:,1);
                    j=j+1;
                    clear planes planesTemp
                    for k = 1:raw3D.l
                        if ismember(k,searched)==1
                            planesTemp(:,1) = obj.nbors(k,1:nnz(obj.nbors(k,:)));
                            planes = planesTemp;
                            searched(searched==k) = [];
                            break
                        end
                        if k == raw3D.l
                            working = 0;
                        end
                    end
                end
            end
            obj.raw = obj.raw;
            
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