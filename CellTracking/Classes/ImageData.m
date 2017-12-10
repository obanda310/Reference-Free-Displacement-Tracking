classdef ImageData
    properties
        Black
        Trans
        Fluor
        Area
        Borders
        ROIstack
    end
    methods
        function obj = ImageData(autoChk)
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
            
            %Open processed image stack of dots
            files = dir('*.tif'); %Check Directory for default filenames
            for k = 1:length(files)
                current=files(k).name;
                check(k)=strcmp(current(end-6:end),'roi.tif');
            end
            loc=find(check);
            if size(loc,1)==1
                obj.ROIstack = getImages(files(loc(1)).name);
            else
                obj.ROIstack = getImages();
            end
            
            obj.Borders = ones(size(obj.Area,1),size(obj.Area,2));
        end
    end
end