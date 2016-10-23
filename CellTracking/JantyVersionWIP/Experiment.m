classdef Experiment
    properties
        images
        metadata
        cellImg
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
        end
        function [roiStack,roiCell,roiBounds,bkImg] = cropImgs(obj)
            roiGUIHandle = ROI_GUI(obj);
            waitfor(roiGUIHandle);
            
            roiStack = roiGUIData.ROI;
            roiBounds = roiGUIData.ROIBounds;
            
            roiCell = imcrop(exp.cellImg,roiBounds);
            
            bkImg = zeroes(roiBounds(4)+1,roiBounds(3)+1);
        end
    end
end