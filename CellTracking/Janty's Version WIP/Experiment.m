classdef Experiment
    properties
        trajImg     % image of the trajectories map output from imageJ
        cellImg     % image of the cell(s) being observed
        trajData    % matrix of trajectory data
        numFrames   % number of frames (z-slices) in experiment
        numTrajs    % number of trajectories in experiment
        resolution  % [y resolution, x resolution]
    end
    methods
        function obj = Experiment(trajImgFile,cellImgFile,trajDataFile)
            if nargin == 0
                [tFile,tPath] = uigetfile('*.tif','Select trajectories image');
                trajImgFile = [tPath,tFile];
                [cFile,cPath] = uigetfile('*.tif','Select cell overlay image');
                cellImgFile = [cPath,cFile];
                [dFile,dPath] = uigetfile('*.xlsx');
                trajDataFile = [dPath,dFile];
            end
            obj.trajImg = imread(trajImgFile);
            obj.cellImg = imread(cellImgFile);
            [obj.trajData,~,~] = xlsread(trajDataFile);
            
            obj.resolution = size(obj.trajImg);
            % Subtract y-position data from y-resolution to redefine the
            % origin in the top-left corner of the image instead of the
            % bottom-left corner.
            obj.trajData(:,5) = obj.resolution(2)-obj.trajData(:,5);
            
            obj.numFrames = max(obj.trajData(:,3)) + 1;
            obj.numTrajs = max(obj.trajData(:,2));
        end
    end
end