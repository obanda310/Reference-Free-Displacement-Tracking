% Open czi file and extract relevant metadata
function [images,meta] = getCzi(file)
    %% Import image stack
    if nargin == 1
        filename = file;
    elseif nargin == 0
        % Choose .czi image stack
        [name,path] = uigetfile('*.czi');
        filename = [path,name];
    end
    
    % Add the functions in the bfmatlab folder and its subfolders to the
    % path available for this code
    addpath(genpath('bfmatlab'));
    
    % Read into matlab
    czi = bfopen(filename);
    % Get image dimensions from first image in .czi file
    [imgRows, imgCols] = size(czi{1}{1,1});
    % Images are stored in the first cell in the czi cell array
    imageData = czi{1,1};
    % Get number of images (z-slices) from .czi file 
    noImgs = length(imageData);
    % Preallocate 3D image matrix
    images = zeros(imgRows,imgCols,noImgs);
    % Store images in one 3D image matrix (instead of czi cell array)
    parfor i = 1:noImgs
        % Column 1 in imageData contains the image grayscale intensity 
        % matrices; row index i corresponds to position in z-stack
        images(:,:,i) = imageData{i,1};
    end
    %% Read metadata
    % Metadata are stored in the second cell in the czi cell array
    metadata = czi{1,2};   
    % Need to use ".get" command to access metadata since the data is
    % stored in a hashtable. Ctrl+F "matlab" in the reference:
    % "bfmatlab_metadata_reference.pdf" for more information
    % To see all available metadata easily, see "all_metadata.xlsx"
    meta.binning          = metadata.get('Global Information|Image|Channel|Binning #1');
    meta.zStart           = str2double(metadata.get('Global Experiment|AcquisitionBlock|MultiTrackSetup|ZStackSetup|First|Distance|Value #1')); %um
    meta.zEnd             = str2double(metadata.get('Global Experiment|AcquisitionBlock|MultiTrackSetup|ZStackSetup|Last|Distance|Value #1')); %um
    meta.exposureTime     = metadata.get('Global HardwareSetting|ParameterCollection|ExposureTime #1');
    meta.colorDepth       = str2double(metadata.get('Global Information|Image|ComponentBitCount #1'));
    meta.scaling          = str2double(metadata.get('Global Scaling|Distance|Value #1'));
    meta.sizeX            = str2double(metadata.get('Global Information|Image|SizeX #1'));
    meta.sizeY            = str2double(metadata.get('Global Information|Image|SizeY #1'));
end