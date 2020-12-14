classdef RawData2D %replacing Shear data
    properties
        pxRaw
        raw
        col
        Int
        rawX %1 = Raw X Centroid Location of Traj in First Frame
        rawY %2 = Raw Y Centroid Location of Traj in First Frame
        rawX1
        rawY1
        rawZ
        firstFrame %3 = First Frame that a Traj Appears
        lastFrame %4 = Last Frame that a Traj Appears
        lastX
        lastY
        numTraj
        rawdX
        rawdY
        rawdXY
        ID
    end
    methods
        function [obj] = RawData2D(image,raw)
            
            % Sub-pixel Object Detection Using Kilfoil 'feature2D' Function
            clear xyz2D
            %[r] = feature2D(img,lambda,w,masscut,Imin)
            for i = 1:size(image.MaskStack,3)
                clear temp
                currentImg = image.MaskStack(:,:,i);
                temp = feature2D(currentImg,1,3,2000,50);
                if i == 1
                    temp(:,6) = i;
                    xyz2D = temp;
                else
                    temp(:,6) = i;
                    xyz2D = cat(1,xyz2D,temp);
                end
            end
            
            % Here we will link objects in the Z direction
            disp('Linking 3D Detections Between Planes.')
            % Set a maximum linking distance in microns that any object can still be
            % considered part of a pillar. Smaller values will speed up code.
            maxLinkDistance = 1.2; % !!**value should be smaller than array spacing**!!
            maxLD = maxLinkDistance/raw.dataKey(9,1);
            disp(['Maximum allowable XY Translation (Microns): ',num2str(maxLinkDistance)])
            % Set a maximum number of frames to look for a linked object before giving
            % up (maxJumpDistance)
            maxJD = 15;
            % Linking using trackmem function (kilfoil code)
            disp('Linking Pillars')
            [lub] = trackmem(xyz2D,maxLD,2,0,maxJD);
            obj.pxRaw = lub;
            lub(:,1:2) = lub(:,1:2)*raw.dataKey(9,1);
            lub(:,6) = lub(:,6)*raw.dataKey(10,1);
            obj.raw = lub;
            
            
            
            %%Object Detection Plot
            detections = figure;
            imshow(image.Black)
            hold on
            map = brewermap(max(obj.pxRaw(:,6)),'*Spectral');
            
            for i = 1:max(obj.pxRaw(:,6)) %size(subpixMaxima,3)-10
                [tempFrame,~] = find(obj.pxRaw(:,6) == i);
                scatter3(obj.pxRaw(tempFrame,1),obj.pxRaw(tempFrame,2),obj.pxRaw(tempFrame,6),'.','SizeData',100,'MarkerFaceColor',map(i,1:3),'MarkerEdgeColor',map(i,1:3))
            end
            hold off
            filePath = cd;
            savefile = [filePath '\Tracking_Unlinked Detections.tif'];
            export_fig(detections,savefile,'-native');
        end
        function obj = shape2DData(obj,raw)
            disp('Copying Raw data Into ShearData Class Variable')
            progressbar('Indexing displacement data')
            
            obj.numTraj = max(obj.raw(:,7)); %Number of Trajectories
            %HANDLING XYZ DATA
            skipCount = 0;
            missingNo = 0;
            for i = 1:obj.numTraj
                progressbar(i/obj.numTraj)
                % Here we build a book of pages (3D array) with the data for a single
                % object/trajectory per page. Most of the code is to ensure that each
                % page is the same size matrix as the next.
                
                % stores data to 'tempObj' pertaining to all frames of one
                % object/trajectory
                tempObj = obj.pxRaw(obj.pxRaw(:,7)==i,:);
                % number of frames that the current object appears in.
                numFrames = size(tempObj,1);
                
                % the first frame that an object appears in (tracking software starts
                % at 0
                obj.firstFrame(i) = min(tempObj(:,6));
                % the last frame that an object appears in
                obj.lastFrame(i) = max(tempObj(:,6));
                for j = 1:size(tempObj,1)
                obj.rawX(tempObj(j,6),i) = tempObj(j,1);
                obj.rawY(tempObj(j,6),i) = tempObj(j,2);
                obj.rawZ(tempObj(j,6),i) = tempObj(j,6);
                obj.Int(tempObj(j,6),i) = tempObj(j,3);
                end
                
            end
            obj.numTraj = size(obj.rawX,2);
            
            disp('Storing XY Displacements')
            for i = 1:obj.numTraj
                obj.rawX1(i) = obj.rawX(obj.firstFrame(i),i); %initial x
                obj.rawY1(i) = obj.rawY(obj.firstFrame(i),i); %initial y
                obj.rawdX(obj.firstFrame(i):obj.lastFrame(i),i) = obj.rawX(obj.firstFrame(i):obj.lastFrame(i),i) - obj.rawX1(i); %frame specific dx
                obj.rawdY(obj.firstFrame(i):obj.lastFrame(i),i) = obj.rawY(obj.firstFrame(i):obj.lastFrame(i),i) - obj.rawY1(i); %frame specific dy
                obj.rawdXY(obj.firstFrame(i):obj.lastFrame(i),i) = (obj.rawdX(obj.firstFrame(i):obj.lastFrame(i),i).^2 + obj.rawdY(obj.firstFrame(i):obj.lastFrame(i),i).^2).^0.5; %frame specific total displacement
                obj.ID(i) = i; %pillar ID
            end
            
            
        end
        %%
        function obj = link3D(obj,r,raw) % ***currently not in use***
            %=% Link 3D detections to 2D Detections and Apply Corrections
            % 2D Centroids (in XY) can provide locations of 'incomplete' features near the
            % surface, where deformations are highest, and where 3D object detection can fail
            % to detect objects. Unfortunately, 2D info is subject to additional flaws,
            % such as tilted ellipsoids (caused by imperfect laser alignment), or
            % tilted columns (caused by misaligned patterning plane, imaging plane,
            % and hydrogel surface. This section will attempt to use 3D info as a
            % guide to catalogue 2D info.
            
            xyz2Ds = obj.raw;
            xyz2Ds(:,3:5) = [];
            xyz2Ds(:,1) = xyz2D(:,2);
            xyz2Ds(:,2) = xyz2D(:,1);
            xyz2Ds(:,1:2) = xyz2Ds(:,1:2)*raw.dataKey(9,1); %convert to microns
            xyz2Ds(:,3) = xyz2Ds(:,3)*raw.dataKey(10,1); %convert to microns
            % Preserve original 2D list
            xyz2Di = xyz2Ds;
            xyz2Ds = xyz2Di; %use to reset
            % Creat a final list with 3D assignments
            xyz2Dfinal = xyz2Ds;
            
            % First link 3D objects to the 2D objects within a maximum allowable
            % distance
            maxLD = 1.5; %should be less than half of the Z spacing during patterning
            clear xyz2Dci xyz2Ddi rDe
            rDe = delaunayn(double([r.X, r.Y, r.Z]));
            [xyz2Dci,xyz2Ddi] = dsearchn(double([r.X, r.Y, r.Z]),rDe,xyz2Ds); %indices and distance of closest 2D detections
            xyz2Dci(xyz2Ddi(:,1)>maxLD,1) = 0;
            xyz2Ddi(xyz2Ddi(:,1)>maxLD,1) = 0;
            xyz2Dfinal(:,4) = xyz2Dci(:,1);
            xyz2Dfinal(:,5) = xyz2Ddi(:,1);
            xyz2Dfinal(xyz2Dfinal(:,4)>0,6) = r.Z(xyz2Dfinal(xyz2Dfinal(:,4)>0,4))-xyz2Dfinal(xyz2Dfinal(:,4)>0,3);
            xyz2Dfinal = sortrows(xyz2Dfinal,[4 6]);
            
            lowCO = raw.dataKey(10,1) - maxLD;
            highCO = maxLD - raw.dataKey(10,1);
            maxLD = 1;
            go =1;
            while go == 1
                clear pEnds
                sizeBef = size(xyz2Ds,1);
                xyz2Ds = xyz2Dfinal(xyz2Dfinal(:,4)==0,1:3);
                sizeAft = size(xyz2Ds,1);
                if sizeBef == sizeAft
                    break
                end
                pEnds = cat(1,xyz2Dfinal(xyz2Dfinal(:,6)<lowCO,1:3),xyz2Dfinal(xyz2Dfinal(:,6)>highCO,1:3));
                pIdx = cat(1,xyz2Dfinal(xyz2Dfinal(:,6)<lowCO,4),xyz2Dfinal(xyz2Dfinal(:,6)>highCO,4));
                try
                    rDe = delaunayn(pEnds);
                catch
                    break
                end
                [xyz2Dci,xyz2Ddi] = dsearchn(pEnds,rDe,xyz2Ds);
                xyz2Dci(xyz2Ddi(:,1)>maxLD,1) = 0;
                xyz2Ddi(xyz2Ddi(:,1)>maxLD,1) = 0;
                xyz2Didx = zeros(size(xyz2Dci));
                xyz2Didx(xyz2Dci(:,1)>0,1) = pIdx(xyz2Dci(xyz2Dci(:,1)>0),1);
                xyz2Dfinal(1:size(xyz2Dci,1),4) = xyz2Didx(:,1);
                xyz2Dfinal(1:size(xyz2Ddi,1),5) = xyz2Ddi(:,1);
                xyz2Dfinal(xyz2Dfinal(:,4)>0,6) = r.Z(xyz2Dfinal(xyz2Dfinal(:,4)>0,4),1)-xyz2Dfinal(xyz2Dfinal(:,4)>0,3);
                xyz2Dfinal = sortrows(xyz2Dfinal,[4 6]);
                lowCO = lowCO - raw.dataKey(10,1);
                highCO = highCO + raw.dataKey(10,1);
                
            end
            
            figure
            scatter3(xyz2Dfinal(xyz2Dfinal(:,4)==0,1),xyz2Dfinal(xyz2Dfinal(:,4)==0,2),xyz2Dfinal(xyz2Dfinal(:,4)==0,3))
            % Then use the farthest assigned 2D detections to continue assignments
            % until no more assignments are possible
            
            % Iteratively assign a 3D object index to a 2D object index if it is
            % within the allowable distance, and remove that 2D object from the search
            % list. Then repeat the search
           
            
            
            %--------------------------------------------------------
            % Calculate average slope of pillars in Z
            
%             for i = 1:size(r.X,1)
%                 if r.colSize(i) >1 && ismember(i,r.ND)
%                     clear tempcol
%                     tempcol(:,1:3) = r.r(r.col(:)==r.col(i),1:3);
%                     sortrows(tempcol,3);
%                     Xs(i) = tempcol(1,1)-tempcol(end,1);
%                     Ys(i) = tempcol(1,2)-tempcol(end,2);
%                     Zs(i) = tempcol(1,3)-tempcol(end,3);
%                 end
%             end
%             colV(1) = mean(Xs);
%             colV(2) = mean(Ys);
%             colV(3) = mean(Zs);
%             colV = colV/min(colV);
        end
    end
end