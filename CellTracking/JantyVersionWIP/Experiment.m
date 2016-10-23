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
            % If a cellFile input argument is provided, use that as an
            % input for getImages without asking the user to select a file
            % manually
            else
                [obj.cellImg,~] = getImages(cellFile);
            end
            
            % Process images and get centroid locations
            [masks,obj.centroids2d] = ... 
                getCentroidsStack(obj.images,obj.metadata);
            % Clear any objects that are not 3-Dimensionally larger than a 
            % threshold value of 20 pixels
            obj.masks = bwareaopen(masks,20); 
        end
        function [roiImgs,roiMasks,roiCell,roiBounds,bkImg] = cropImgs(obj)
            roiGUIHandle = ROI_GUI(obj);
            w = msgbox('Click ok when you''re happy with the ROI');
            waitfor(w);
            roiGUIData = getappdata(roiGUIHandle);
            close(roiGUIHandle);

            roiImgs = roiGUIData.ROI;
            roiBounds = roiGUIData.ROIBounds;
            
            noImgs = size(obj.images,3);
            roiMasks = zeros(roiBounds(4)+1,roiBounds(3)+1,noImgs);
            for i = 1:noImgs
                newImg = imcrop(obj.masks(:,:,i),roiBounds);
                roiMasks(:,:,i) = newImg;
            end
            
            roiCell = imcrop(obj.cellImg,roiBounds);
            
            bkImg = zeros(roiBounds(4)+1,roiBounds(3)+1);
        end
    end
end