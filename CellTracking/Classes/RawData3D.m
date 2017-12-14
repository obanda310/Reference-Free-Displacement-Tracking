classdef RawData3D
    properties
        s
        r
        X
        Y
        Z
        l
        Mass
        MaxI
        rg
        PctAbove
        row %8
        col %9
        XSC %10 Shift correction X
        YSC %11 Shift Correction Y
        ND %non-deformed markers
        
    end
    methods
        
        function obj = RawData3D(res,raw,r)
            if nargin == 2
                obj.s(1,1) = size(res,1); obj.s(2,1) = size(res,2); obj.s(3,1) = size(res,3);
                obj.s(1:2,2) = obj.s(1:2,1)*raw.dataKey(9,1);
                obj.s(3,2) = obj.s(3,1)*raw.dataKey(10,1);
                d = round(1/raw.dataKey(9,1));
                dz = round(1.5/raw.dataKey(10,1));
                obj.r=...
                    feature3dMB(res, d , [d d dz], [obj.s(1) obj.s(2) obj.s(3)],[1 1 1],round(.4/raw.dataKey(9,1)),0,.05); %
                obj.r(:,1:2) = obj.r(:,1:2)*raw.dataKey(9,1);
                obj.r(:,3) = obj.r(:,3)*raw.dataKey(10,1);
                obj.l = size(obj.r,1);
            elseif nargin == 3
                obj.s(1,1) = size(res,1); obj.s(2,1) = size(res,2); obj.s(3,1) = size(res,3);
                obj.s(1:2,2) = obj.s(1:2,1)*raw.dataKey(9,1);
                obj.r = r;
                obj.l = size(r,1);
            end
        end
        
        function obj = TranscribeR(obj)
            obj.X = obj.r(:,1);
            obj.Y = obj.r(:,2);
            obj.Z = obj.r(:,3);
            obj.Mass = obj.r(:,4);
            obj.rg = obj.r(:,5);
            obj.MaxI = obj.r(:,6);
            obj.PctAbove = obj.r(:,7);
        end
        
        function viewDetections(obj,raw,offset,image)
            figure
            if nargin == 2
            scatter3(obj.X,obj.s(2)*raw.dataKey(9,1)-obj.Y,obj.Z)
            elseif nargin == 4
                imshow(image)
                hold on
                scatter3((obj.X-offset(1,1))/raw.dataKey(9,1),((obj.Y)-offset(1,2))/raw.dataKey(9,1),obj.Z)
                hold off
            end
        end
        
        function obj = regionCheck(obj,image,raw)
            obj.ND = zeros(1,1);
            for i = 1:obj.l
                %if it is under the cell
                if image.ADil(round(obj.Y(i)/raw.dataKey(9,1)),round(obj.X(i)/raw.dataKey(9,1)))~=0
                    obj.ND = cat(1,obj.ND,i);
                end
            end
            obj.ND(1,:) = [];
        end
    end
end
