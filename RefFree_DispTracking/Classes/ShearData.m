classdef ShearData
    properties
        
        %Index List for Book1 - Frame Dependent Values
        rawX %1 = Raw X Centroid Location of Traj in Specified Frame
        rawY %2 = Raw Y Centroid Location of Traj in Specified Frame
        rawZ % Frame position in microns
        rawdX %3 = Raw dX Centroid Location of Traj in Specified Frame from First Frame
        rawdY %4 = Raw dY Centroid Location of Traj in Specified Frame from First Frame
        rawdXY %5 = Magnitude of displacement
        Int %6 = Intensity in current frame
        GInt %7 = Gauss filtered intensity values (created/used in section 4)
        gtdX %8 = Global Tilt Corrected dX
        gtdY %9 = Global Tilt Corrected dY
        gtdXY %10 = Global Tilt Corrected Magnitude
        mltdX %11 = Mean local dx tilt
        mltdY %12 = Mean local dy tilt
        ltdX %13 = Local Tilt Corrected dx
        ltdY %14 = Local Tilt Corrected dy
        ltdXY %15 = Local Tilt Corrected Magnitude
        coCheck %16 = Index %15 is above cmCutoff (section 4.6)
        coltdX %17 = Local Tilt Corrected dx above cutoff (using idx 16 as a mask)
        coltdY %18 = Local Tilt Corrected dy above cutoff
        coltdXY %19 = Local Tilt Corrected Magnitude above cutoff
        
        %Index List for Book2 - Frame Independent Values
        rawX1 %1 = Raw X Centroid Location of Traj in First Frame
        rawY1 %2 = Raw Y Centroid Location of Traj in First Frame
        firstFrame %3 = First Frame that a Traj Appears
        lastFrame %4 = Last Frame that a Traj Appears
        lastdX%5 = Value of book1 index 3 (see above) in Last Frame that Traj Appears
        lastdY%6 = Value of book1 index 4 (see above) in Last Frame that Traj Appears
        lastX%7 = Value of book1 index 1 (see above) in Last Frame that Traj Appears
        lastY%8 = Value of book1 index 2 (see above) in Last Frame that Traj Appears
        MaxdXY%9 = Maximum Magnitude of displacement
        ID %10= Numeric ID of Pillar
        gtLastdX%11= Global Tilt Corrected Final dx
        gtLastdY%12= Global Tilt Corrected Final dy
        gtLastdXY%13= Global Tilt Corrected Final Magnitude
        Top1%14= Predicted Top Surface
        dTop1%15= Deviation from top Surface
        Top2%16= Alternative Predicted Top Surface
        dTop2%17= Deviation from Alternative Top Surface
        ltLastdX%18= Local Tilt Corrected Final dx
        ltLastdY%19= Local Tilt Corrected Final dy
        ltLastdXY%20= Local Tilt Corrected Final Magnitude
        coCheck2%21= Book1 Index 16 is >0 for at least 3 consecutive frames
        coFilt %22
        
        %Other Variables
        numTraj % total number of columns
        numFrames % total number of frames
        noCellTraj
        noiseBook %
        cutoff
        cmCutoff
        sumIndFinal
        
    end
    methods
        %Currently:
        %ObjectShearData - Loads raw tracking data and calculates displacement
        %information
        
        %globalTilt - Adds a correction factor to displacement data based
        %on an approximate global tilt of the sample
        
        %localTilt - Refines tilt measures using local appromiximations of
        %tilt
        
        %lTcutoff - Thresholds data using a user input cutoff (fast way to
        %remove noisy regions).
        
        %% ------------------------------------------------------------
        function obj = ShearData(raw,image)
            disp('Copying Raw data Into ShearData Class Variable')
            progressbar('Indexing displacement data')
            
            obj.numTraj = max(raw.data(:,raw.dataKey(4,1))); %Number of Trajectories
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
                tempObj = raw.data(raw.data(:,raw.dataKey(4,1))==i,:);
                % number of frames that the current object appears in.
                numFrames = size(tempObj,1);
                if numFrames > 0
                    % the first frame that an object appears in (tracking software starts
                    % at 0
                    startFrame = min(tempObj(:,raw.dataKey(3,1)))+raw.dataKey(8,1);
                    % the last frame that an object appears in
                    endFrame = (startFrame+numFrames-1);
                    
%                     if image.Borders(ceil(tempObj(1,3)*raw.dataKey(7,1)),ceil(tempObj(1,2)*raw.dataKey(7,1))+1)==0
%                         skip = 1;
%                         skipCount = skipCount+1
%                     else
                        skip = 0;
%                     end
                    
                    if skip == 0
                        obj.firstFrame((i-skipCount)-missingNo) = startFrame;
                        obj.lastFrame((i-skipCount)-missingNo) = endFrame;
                        % this fills in the upper portion of obj matrix with zeros if the first
                        % frame is not 0
                        
                        obj.rawX(startFrame:endFrame,(i-skipCount)-missingNo) = tempObj(:,raw.dataKey(1,1)).*raw.dataKey(7,1);
                        obj.rawY(startFrame:endFrame,(i-skipCount)-missingNo) = tempObj(:,raw.dataKey(2,1)).*raw.dataKey(7,1);
                        obj.rawZ(startFrame:endFrame,(i-skipCount)-missingNo) = tempObj(:,9).*raw.dataKey(7,1);
                        obj.Int(startFrame:endFrame,(i-skipCount)-missingNo) = tempObj(:,raw.dataKey(5,1));
                    end
                else
                    missingNo = missingNo +1;
                end
                if i == 1 || i == 5000 || i == 10000 || i == 15000 || i==20000 || i==25000
                    disp(['Progress: ' num2str(i) ' of ' num2str(obj.numTraj)])
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
        %%-------------------------------------------------------------
        
        
        
        
        
        
        %% ------------------------------------------------------------
        function [obj] = globalTilt(obj,image,raw,filePath)
            obj.numFrames = size(obj.rawX,1);
            % Find non-deformed pillars with no previous info
            noCellTrajIni = 0;
            SE = strel('disk',80);
            imageAreaDil = imerode(image.Area,SE);
            
            for i = 1:obj.numTraj
                
                %if it is in black region
                if imageAreaDil(ceil(obj.rawY1(i)/raw.dataKey(9,1)),ceil(obj.rawX1(i)/raw.dataKey(9,1)))~=0
                    if noCellTrajIni == 0
                        noCellTrajIni = i;
                    else
                        noCellTrajIni = cat(1,noCellTrajIni,i);
                    end
                end
            end
            figure
            imshow(imageAreaDil)
            
            % Tilt Correction
            [obj.noiseBook,noiseStats,obj.sumIndFinal] = tiltCorrection(image.ROIstack,image.Trans,obj,noCellTrajIni);
            
            
            %
            for i = 1:obj.numFrames
                obj.gtdX(i,:) = ((obj.rawdX(i,:)) - obj.noiseBook(i,2)).*(obj.rawdX(i,:)~=0);
                obj.gtdY(i,:) = ((obj.rawdY(i,:)) - obj.noiseBook(i,3)).*(obj.rawdY(i,:)~=0);
                obj.gtdXY(i,:) = ((obj.gtdX(i,:).^2)+(obj.gtdY(i,:).^2)).^.5;
                obj.noiseBook(i,6) = mean(obj.gtdX(i,obj.sumIndFinal));
                obj.noiseBook(i,7) = mean(obj.gtdY(i,obj.sumIndFinal));
                obj.noiseBook(i,8) = mean(obj.gtdXY(i,obj.sumIndFinal));
            end
            
            for i = 1:obj.numTraj
                
                obj.lastdX(i) = obj.rawdX(obj.lastFrame(i),i);
                obj.lastdY(i) = obj.rawdY(obj.lastFrame(i),i);
                obj.lastX(i) = obj.rawX(obj.lastFrame(i),i);
                obj.lastY(i) = obj.rawY(obj.lastFrame(i),i);
                obj.MaxdXY(i) = max(obj.rawdXY(:,i));%(obj.lastdX(i).^2 + obj.lastdY(i).^2).^0.5;
                obj.gtLastdX(i) = obj.gtdX(obj.lastFrame(i),i);
                obj.gtLastdY(i) = obj.gtdY(obj.lastFrame(i),i);
                obj.gtLastdXY(i) = (obj.gtLastdX(i).^2 + obj.gtLastdY(i).^2).^0.5;
            end
            
            % Find all pillars in non-deformed regions
            clear imageArea2
            imageArea2 = zeros(size(image.Area,1)+100,size(image.Area,2)+100);
            imageArea2(51:size(image.Area,1)+50,51:size(image.Area,2)+50) = image.Area;
            
            %Add areas of high deformation to the binary mask of cell to
            %ignore these areas for tilt correction
%             for i = 1:obj.numTraj
%                 if obj.gtLastdXY(i) > (1)
%                     imageArea2((round(obj.rawY1(i)/raw.dataKey(9,1))):(round(obj.rawY1(i)/raw.dataKey(9,1))+100),(round(obj.rawX1(i)/raw.dataKey(9,1))):(round(obj.rawX1(i)/raw.dataKey(9,1))+100))=0;
%                 end
%             end
            
            imageArea2 = imcrop(imageArea2,[51,51,size(image.Area,2)-1,size(image.Area,1)-1]);
            obj.noCellTraj = 0;
            for i = 1:obj.numTraj
                %if it is in black region
                if imageArea2(round(obj.rawY1(i)/raw.dataKey(9,1)),round(obj.rawX1(i)/raw.dataKey(9,1)))~=0
                    if obj.noCellTraj == 0
                        obj.noCellTraj = i;
                    else
                        obj.noCellTraj = cat(1,obj.noCellTraj,i);
                    end
                end
            end
            figure
            imshow(imageArea2)
            
            %CHANGE FONT SIZES HERE
            AxisFontSize = 28;
            AxisTitleFontSize = 28;
            LegendFontSize = 20;
            %
            mkdir('Histograms')
            bins = 0:.025:3;
            errorHist = figure;
            hold on
            disp('obj.noCellTraj')
            obj.noCellTraj
            histogram(obj.gtdXY(:,obj.noCellTraj),bins,'normalization','probability')
            histogram(obj.rawdXY(:,obj.noCellTraj),bins,'normalization','probability')
            set(gca,'fontsize',AxisFontSize)
            xt = 'Reference Error';% input('enter the xaxis label','s');
            yt = 'Probability'; %input('enter the yaxis label','s');
            xl = xlabel(xt);
            yl = ylabel(yt);
            le2 = 'Normal'; %input('enter the legend','s');
            le = 'Global Tilt';
            leg = legend(le,le2,'location','northeast');
            leg.FontSize = LegendFontSize;
            axis([0 1 0 .4])
            savefile = [filePath '\Histograms' '\GTErrorHistNoStress.tif'];
            export_fig(errorHist,savefile);
            
            
            errorHist2 = figure;
            hold on
            histogram((obj.gtdXY(:,:)),bins,'normalization','probability')
            histogram((obj.rawdXY(:,:)),bins,'normalization','probability')
            set(gca,'fontsize',AxisFontSize)
            xt = 'Reference Error';% input('enter the xaxis label','s');
            yt = 'Probability'; %input('enter the yaxis label','s');
            xl = xlabel(xt);
            yl = ylabel(yt);
            le2 = 'Normal'; %input('enter the legend','s');
            le = 'Global Tilt';
            leg = legend(le,le2,'location','northeast');
            leg.FontSize = LegendFontSize;
            axis([0 1 0 .4])
            savefile = [filePath '\Histograms' '\GTErrorHistAll.tif'];
            export_fig(errorHist2,savefile);
            %
            
            
        end
        
        %%-------------------------------------------------------------
        
        
        
        
        
        
        %% ------------------------------------------------------------
        function obj = localTilt(obj,filePath)
            
            clear book4
            %book4 contains the closest 500 pillars to index pillar(dim 1) that do not
            %fall within the boundaries of the cell specified in var 'image.Area'
            
            % find mean deviation of 500 closest in 'noCellTraj' variable
            for i = 1:obj.numTraj
                clear tempDistances
                if size(obj.noCellTraj,1)>100
                    tempDistances(1:size(obj.noCellTraj,1),1) = ((obj.rawX1(obj.noCellTraj)-obj.rawX1(i)).^2.+(obj.rawY1(obj.noCellTraj)-obj.rawY1(i)).^2).^0.5;
                    [tempDistances2,origOrder] = sortrows(tempDistances);
                    book4(i,1:100) = obj.noCellTraj(origOrder(1:100,1),1); %record the closest 500 pillars
                else
                    tempDistances(1:size(obj.noCellTraj,1),1) = ((obj.rawX1(:)-obj.rawX1(i)).^2.+(obj.rawY1(:)-obj.rawY1(i)).^2).^0.5;
                    [~,origOrder] = sortrows(tempDistances);
                    book4(i,1:size(obj.noCellTraj,1)) = origOrder(1:size(obj.noCellTraj,1),1); %record all pillars
                end
            end
            clear tempDistances
            
            tempX = obj.rawdX;
            tempY = obj.rawdY;
            tempX(tempX==0) = NaN;
            tempY(tempY==0) = NaN;
            
            for j = 1:obj.numFrames
                for i = 1:obj.numTraj
                    if obj.rawdX(j,i) == 0 || obj.rawdY(j,i) == 0
                        obj.mltdX(j,i) = 0;
                        obj.mltdY(j,i) = 0;
                        obj.ltdX(j,i) = 0;
                        obj.ltdY(j,i) = 0;
                        obj.ltdXY(j,i) = 0;
                    else
                        obj.mltdX(j,i) = nanmean(tempX(j,book4(i,:)));
                        obj.mltdY(j,i) = nanmean(tempY(j,book4(i,:)));
                        obj.ltdX(j,i) = obj.rawdX(j,i)-obj.mltdX(j,i);
                        obj.ltdY(j,i) = obj.rawdY(j,i)-obj.mltdY(j,i);
                        obj.ltdXY(j,i) = ((obj.ltdX(j,i).^2)+(obj.ltdY(j,i).^2)).^.5;
                    end
                end
            end
            
            
            
            for i = 1:obj.numTraj
                
                obj.ltdX(isnan(obj.ltdX)) = 0;
                obj.ltdY(isnan(obj.ltdY)) = 0;
                obj.ltLastdX(i) = obj.ltdX(find(abs(obj.ltdX) == max(abs(obj.ltdX(:,i))),1,'last'));%obj.ltdX(obj.lastFrame(i),i);
                obj.ltLastdY(i) = obj.ltdY(find(abs(obj.ltdY) == max(abs(obj.ltdY(:,i))),1,'last'));%obj.ltdY(obj.lastFrame(i),i);
                obj.ltLastdXY(i) = max(obj.ltdXY(:,i));%(obj.ltLastdX(i).^2 + obj.ltLastdY(i).^2).^0.5;
                
            end
            % for i = 1:obj.numTraj
            %     if obj.ltdX(obj.lastFrame(i)-1,i) == 0 || obj.ltdY(obj.lastFrame(i)-1,i) == 0
            %     obj.ltLastdX(i) = max(obj.ltdX(:,i));
            %     obj.ltLastdY(i) = max(obj.ltdY(:,i));
            %     obj.ltLastdXY(i) = (obj.ltLastdX(i).^2 + obj.ltLastdY(i).^2).^0.5;
            %     else
            %     obj.ltLastdX(i) = obj.ltdX(obj.lastFrame(i)-1,i);
            %     obj.ltLastdY(i) = obj.ltdY(obj.lastFrame(i)-1,i);
            %     obj.ltLastdXY(i) = (obj.ltLastdX(i).^2 + obj.ltLastdY(i).^2).^0.5;
            %     end
            % end
            
            
            %CHANGE FONT SIZES HERE
            AxisFontSize = 28;
            AxisTitleFontSize = 28;
            LegendFontSize = 20;
            %
            mkdir('Histograms')
            bins = 0:.025:3;
            errorHist = figure;
            hold on
            histogram(obj.ltdXY(:,obj.noCellTraj),bins,'normalization','probability')
            histogram(obj.gtdXY(:,obj.noCellTraj),bins,'normalization','probability')
            histogram(obj.rawdXY(:,obj.noCellTraj),bins,'normalization','probability')
            set(gca,'fontsize',AxisFontSize)
            xt = 'Reference Error';% input('enter the xaxis label','s');
            yt = 'Probability'; %input('enter the yaxis label','s');
            xl = xlabel(xt);
            yl = ylabel(yt);
            le = 'Local Tilt';
            le2 = 'Global Tilt';
            le3 = 'Normal'; %input('enter the legend','s');
            leg = legend(le,le2,le3,'location','northeast');
            leg.FontSize = LegendFontSize;
            axis([0 1 0 .4])
            savefile = [filePath '\Histograms' '\ErrorHistNoStress.tif'];
            export_fig(errorHist,savefile);
            
            
            errorHist2 = figure;
            hold on
            histogram((obj.ltdXY(:,:)),bins,'normalization','probability')
            histogram((obj.gtdXY(:,:)),bins,'normalization','probability')
            histogram((obj.rawdXY(:,:)),bins,'normalization','probability')
            set(gca,'fontsize',AxisFontSize)
            xt = 'Reference Error';% input('enter the xaxis label','s');
            yt = 'Probability'; %input('enter the yaxis label','s');
            xl = xlabel(xt);
            yl = ylabel(yt);
            le = 'Local Tilt';
            le2 = 'Global Tilt';
            le3 = 'Normal'; %input('enter the legend','s');
            leg = legend(le,le2,le3,'location','northeast');
            leg.FontSize = LegendFontSize;
            axis([0 1 0 .4])
            savefile = [filePath '\Histograms' '\ErrorHistAll.tif'];
            export_fig(errorHist2,savefile);
        end
        %%-------------------------------------------------------------
        
        
        
        
        
        %% ------------------------------------------------------------
        function obj = lTcutoff(obj)
            obj.cutoff = 0.2;
            obj.cmCutoff = 3; %The first colormap group to be used in filtered data sets
            obj.coCheck(:,:) = obj.ltdXY(:,:)>obj.cutoff;
            for i = 1:obj.numTraj
                if sum(obj.coCheck(:,i))<3
                    obj.coCheck2(i) = 0;
                else
                    obj.coCheck2(i) = 1;
                end
                
                if obj.coCheck2(i) == 1
                    temp = obj.coCheck(:,i);
                    for j = 2:size(temp,1)-1
                        temp2(j,1) = temp(j,1)+temp(j-1,1)+temp(j+1,1);
                    end
                    if max(temp2)<3
                        obj.coCheck2(i) = 0;
                    end
                end
            end
            
            filterSet = find(obj.coCheck2(:)>0);
            obj.coltdX = zeros(obj.numFrames,obj.numTraj);
            obj.coltdY = zeros(obj.numFrames,obj.numTraj);
            obj.coltdXY = zeros(obj.numFrames,obj.numTraj);
            for i = 1:size(filterSet,1);
                obj.coltdX(:,filterSet(i,1)) = obj.ltdX(:,filterSet(i,1));
                obj.coltdY(:,filterSet(i,1)) = obj.ltdY(:,filterSet(i,1));
                obj.coltdXY(:,filterSet(i,1)) = obj.ltdXY(:,filterSet(i,1));
            end
            if size(filterSet,1)<1
                obj.coltdX(:,:) = 0;
                obj.coltdY(:,:) = 0;
                obj.coltdXY(:,:) = 0;
            end
            
            filterMask = obj.coltdXY(:,:)>2;
            obj.coltdX(:,:) = obj.coltdX(:,:).*filterMask;
            obj.coltdY(:,:) = obj.coltdY(:,:).*filterMask;
            obj.coltdXY(:,:) = obj.coltdXY(:,:).*filterMask;
            obj.coFilt(1,:) = obj.coCheck2(:).*obj.ltLastdXY(:);
        end
        %%-------------------------------------------------------------
    end
end