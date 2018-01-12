classdef ImageData
    properties
        Black
        Trans
        Fluor
        Area
        Borders
        ROIstack
        RawStack
        TDots
        ADil
        squares
        SquareBounds
    end
    methods
        function obj = ImageData(autoChk)
            %Load a black image, transmitted image, fluorescent image, and
            %binary image.
            if isa((autoChk),'ImageData')==1
                obj = autoChk;
            else
            obj.Black = loadSpecific('black.tif','*.tif','Select a Black Image of the Correct Dimensions');
            %Open Transmitted Image
            if autoChk == 0
                w = questdlg('Use a Transmitted image for overlays?',...
                    'Transmitted Image (Optional)','Yes','No','Yes');
                waitfor(w);
            else
                w = 'Yes';
            end
            if strcmp(w,'Yes') == 1
                obj.Trans = loadSpecific('Transmitted Cell Image.tif','*.tif','Select Transmitted Image for Overlay');
            else
                obj.Trans = obj.Black;
            end
            
            %Open Fluorescent Image
            if autoChk == 0
                w = questdlg('Use a Fluorescent image for overlays?',...
                    'Fluorescent Image (Optional)','Yes','No','No');
                waitfor(w);
            else
                w = 'No';
            end
            if strcmp(w,'Yes') == 1
                obj.Fluor = loadSpecific('Fluorescent Cell Image.tif','*.tif','Select Fluorescent Image for Overlay');
            else
                obj.Fluor = obj.Black;
            end
            
            %Open Cell Area Binary Image
            if autoChk == 0
                w = questdlg('Input Binary Cell Outline?',...
                    'Binary Outline','Yes','No','Yes');
            else
                w = 'Yes';
            end
            if strcmp(w,'Yes') == 1
                obj.Area = loadSpecific('Binary Mask.tif','*.tif','Select a Thresholded Image of the Cell Area');
                
                if size(obj.Area,1) ~= size(obj.Black,1) || size(obj.Area,2) ~= size(obj.Black,2)
                    [nameAreaFile,filePath] = uigetfile('*.tif','Dimensions of Mask are Incorrect!!!Select a different Thresholded Image of the Cell Area');
                    obj.Area = imread([filePath,nameAreaFile]);
                end
            else
                obj.Area = obj.Black==0;
            end
            
            BT = 35; %pixels 
            obj.Borders = ones(size(obj.Area,1),size(obj.Area,2));
            obj.Borders(1:BT,:) = 0;
            obj.Borders(end-BT:end,:) = 0;
            obj.Borders(:,1:BT) = 0;
            obj.Borders(:,end-BT:end) = 0;
           
                
            end
        end
        
        function obj = roiStack(obj,autoChk)
            %Open processed image stack of markers
            files = dir('*.tif'); %Check Directory for default filenames
            for k = 1:length(files)
                current=files(k).name;
                check(k)=strcmp(current(end-6:end),'roi.tif');
            end
            loc=find(check);
            if size(loc,1)==1 && autoChk ==1
                obj.ROIstack = getImages(files(loc(1)).name);
            else
                obj.ROIstack = getImages();
            end
        end
        
        function obj = rawStack(obj,autoChk)
            %Load the raw fluorecent stack of markers
            matFiles = dir('*.mat');
            if size(matFiles,1)>=1
                for k = 1:length(matFiles)
                    current=matFiles(k).name;
                    if size(current,2)>12
                        check(k)=strcmp(current(end-12:end),'StackName.mat');
                    else
                        check(k) = 0;
                    end
                end
            end
            try
                if size(find(check),2) >0 %&& autoChk ==1
                    load('StackName.mat','StackName')
                    obj.RawStack = single(getImages(StackName));
                else
                    obj.RawStack = single(getImages());
                end
            catch
                obj.RawStack = single(getImages());
            end
        end
        
        function obj = TransDots(obj)
            sumImages = uint16(squeeze(max(permute(obj.RawStack, [3,2,1]))));
            sumImgScale = double(max(max(sumImages)))/(65536);
            sumImages = uint16(sumImages/sumImgScale);
            transImgScale = 65536/mean(prctile(obj.Trans,95));
            obj.Trans = uint16(65536-double((obj.Trans*transImgScale)));
            transImgScale = 65536/mean(prctile(obj.Trans,95));
            obj.Trans = uint16(double((obj.Trans*(transImgScale/2))));
            %imshow(obj.Trans)%invert (should make opaque objects brighter)
            sumImages = sumImages+obj.Trans; %combine dots and cells
            sumImgScale = double(max(max(sumImages)))/65536;
            obj.TDots = uint16(sumImages/sumImgScale);
        end
        
        function obj = DilateBinary(obj,rad)
            SE = strel('disk',rad);
            obj.ADil = imerode(obj.Area,SE); 
        end
        
        function obj = FindNDSquare(obj)
            obj.squares = FindLargestSquares(obj.ADil==255);
            [row,col]=find(obj.squares == max(max(obj.squares)),1,'first');
            length = obj.squares(row,col)-1;
             if length>200
               length = 200;
             end
            obj.SquareBounds(1,1) = col;
            obj.SquareBounds(1,2) = row;
            obj.SquareBounds(1,3) = col+length;
            obj.SquareBounds(1,4) = row+length;
            
        end
    end
end