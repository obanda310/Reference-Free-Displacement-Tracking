% Open czi file and extract relevant metadata
function [images,meta] = getImages(file)
    %% Import image stack
    if nargin == 1
        filename = file;
    elseif nargin == 0
        % Choose image stack
        [name,path] = uigetfile({'*.tif';'*.czi';'*.lsm'});
        filename = [path,name];
    end
    
    % Add the functions in the bfmatlab folder and its subfolders to the
    % path available for this code
    addpath(genpath('bfmatlab'));
    
    % Read into matlab
    file = bfopen(filename);
    % Get image dimensions from first image in file
    [imgRows, imgCols] = size(file{1}{1,1});
    % Images are stored in the first cell in the file cell array
    imageData = file{1,1};
    % Get number of images (z-slices) from file 
    noImgs = size(imageData,1);
    % Preallocate 3D image matrix
    images = zeros(imgRows,imgCols,noImgs);
    % Store images in one 3D image matrix (instead of file cell array)
    for i = 1:noImgs
        % Column 1 in imageData contains the image grayscale intensity 
        % matrices; row index i corresponds to position in z-stack
        images(:,:,i) = imageData{i,1};
    end
    %% Read metadata
    if nargout > 1;
    % Metadata are stored in the second cell in the file cell array
    metadata = file{1,2};
    % Need to use ".get" command to access metadata since the data is
    % stored in a hashtable. Ctrl+F "matlab" in the reference:
    % "bfmatlab_metadata_reference.pdf" for more information
    % To see all available metadata easily, see "all_metadata.xlsx"
    meta.binning          = metadata.get('Global Information|Image|Channel|Binning #1');
    meta.zStart           = str2double(metadata.get('Global Experiment|AcquisitionBlock|MultiTrackSetup|ZStackSetup|First|Distance|Value #1')); %um
    meta.zEnd             = str2double(metadata.get('Global Experiment|AcquisitionBlock|MultiTrackSetup|ZStackSetup|Last|Distance|Value #1')); %um
    meta.exposureTime     = str2double(metadata.get('Global HardwareSetting|ParameterCollection|ExposureTime #1'));
    meta.colorDepth       = str2double(metadata.get('Global Information|Image|ComponentBitCount #1'));
    meta.scalingX         = str2double(metadata.get('Global Scaling|Distance|Value #1'));
    meta.scalingY         = str2double(metadata.get('Global Scaling|Distance|Value #2'));
    meta.scalingZ         = str2double(metadata.get('Global Scaling|Distance|Value #3'));
    meta.sizeX            = size(images,2);
    meta.sizeY            = size(images,1);
    meta.sizeZ            = str2double(metadata.get('Global Information|Image|SizeZ #1'));
    %date                  = metadata.get('Global Information|Image|AcquisitionDateAndTime #1');
    %date = date(1:10);
    meta.date = datetime(date);
    meta.filetype = filename(end-2:end);
    meta.filename = name(1:end-4);
    meta.filepath = path;
    if isnan(meta.scalingX)
        disp('Metadata did not load properly.')
        disp('Please input data manually below as it is required')
        prompt = 'What was the scaling in XY (microns/pixel)? Enter a decimal and press enter: ';
        meta.scalingX = (input(prompt)/1000000);
        meta.scalingY = meta.scalingX;
        prompt = 'What was the color depth of the images? Enter an integer and press enter (ex. 8 for 8bit): ';
        meta.colorDepth = input(prompt); 
    end
    end
    
end