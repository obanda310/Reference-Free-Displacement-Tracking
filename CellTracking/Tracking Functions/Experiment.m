classdef Experiment
    properties
        images
        metadata
        cellImg
        masks
        centroids2d
    end
    methods
        function obj = Experiment(stackFile,cellFile)
            % If a stackFile input argument is not supplied, ask the user
            % to select an image stack he/she would like to analyze
            if ~exist('stackFile','var') || isempty(stackFile)
                w = msgbox('Choose the image stack you''d like to analyze');
                waitfor(w);
                [obj.images,obj.metadata] = getImages;
            % If a stackFile input argument is provided, use that as an
            % input for getImages without asking the user to select a file
            % manually
            else
                [obj.images,obj.metadata] = getImages(stackFile);
            end
            % If a cellFile input argument is not supplied, ask the user
            % to select the cell image he/she would like to use as an
            % overlay
            if ~exist('cellFile','var') || isempty(cellFile)
                w = msgbox('Choose the image of the cells on this gel');
                waitfor(w);
                [obj.cellImg,~] = getImages;
                obj.cellImg = imadjust(uint16(obj.cellImg));
            % If a cellFile input argument is provided, use that as an
            % input for getImages without asking the user to select a file
            % manually
            else
                [obj.cellImg,~] = getImages(cellFile);
                obj.cellImg = imadjust(uint16(obj.cellImg));
            end
            
            % Process images and get centroid locations
            [masks,obj.centroids2d] = ... 
                getCentroidsStack(obj.images,obj.metadata);
            % Clear any objects that are not 2-dimensionally larger than a 
            % threshold value equal to dotSizeThresh
            dotSizeThresh = round(0.9/(1000000*obj.metadata.scalingX)^3);
            obj.masks = bwareaopen(masks,dotSizeThresh); 
        end
        function [roiImgs,roiMasks,roiCell,roiBounds,bkImg] = cropImgs(obj)
            % Use EditStack.m GUI function to select a region of interest
            % and return the bounds of the ROI, as well as the stack of
            % images and masks limited to the selected ROI
            [roiMasks,roiImgs,roiBounds] = ... 
                EditStack(obj.masks,obj.images,obj.centroids2d);
            
            % The ROI image stack is returned as a double, to ensure the
            % images are saved correctly, we convert from double to uint8
            % or uint16, depending on the data type of the original images
            if obj.metadata.colorDepth == 8
                roiImgs = uint8(roiImgs);
            elseif obj.metadata.colorDepth == 16
                roiImgs = uint16(roiImgs);
            end
            
            % Here we save the ROI image stack
            noImgs = size(roiImgs,3);
            % The file path is set such that the ROI image stack will be
            % saved in the same folder as the original image stack, with
            % the same name as the original image stack, but with "roi"
            % appended to the end of the name
            roiFile = [obj.metadata.filepath,obj.metadata.filename,'roi.tif'];
            % If a file already exists with the name we just specified
            % (i.e. if an ROI has been selected before), that file is
            % deleted. This is necessary because the file will not simply
            % be overwritten, because the imwrite function will just append
            % new images to the end of any previously exisiting ones. This
            % is how imwrite handles saving series of images in a single 
            % tiff  file
            if exist(roiFile,'file')
                delete(roiFile)
            end
            % Loop through each image in the ROI image stack, append the 
            % current image (thisImg) to the end of the series, and save it
            % to the filepath specified in roiFile
            for i = 1:noImgs
                % We use imadjust so that the features of the images are
                % visible after being saved
                thisImg = imadjust(roiImgs(:,:,i));
                imwrite(thisImg,roiFile,'WriteMode','append');
            end
            
            % Save the ROI-cropped cell image with a similar file naming
            % convention to that used in saving the ROI stack.
            roiCell = imcrop(obj.cellImg,roiBounds);
            cellFile = [obj.metadata.filepath,obj.metadata.filename,'cell.tif'];
            imwrite(roiCell,cellFile,'tif');
            
            % Save an ROI-sized black image with a similar file naming
            % convention to that used in saving the ROI stack. This image
            % is useful in the trajectories.m function for figure
            % generation purposes
            bkImg = uint16(zeros(roiBounds(4)+1,roiBounds(3)+1));
            blackFile = [obj.metadata.filepath,obj.metadata.filename,'black.tif'];
            imwrite(bkImg,blackFile,'tif');
        end
    end
end