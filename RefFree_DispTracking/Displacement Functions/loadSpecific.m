function [data,filepath] = loadSpecific(nameStr,fileTypeStr,failStr)
%Currently handles only .tif and .txt files

    strL = length(nameStr);
    files = dir(nameStr); %Check Directory for default filenames          
    
    if size(files,1)==1 && strcmp(fileTypeStr,'*.txt')==1
        data = dlmread(files(1).name);
        filepath = files(1).name;
    elseif    size(files,1)==1 && strcmp(fileTypeStr,'*.tif')==1
        data = imread(files(1).name);  
        filepath = files(1).name;
    elseif strcmp(fileTypeStr,'*.txt')==1
        [name,path] = uigetfile(fileTypeStr,failStr);
        file = [path,name];
        data = dlmread(file);
        filepath = file;
    elseif strcmp(fileTypeStr,'*.tif')==1       
        [name,path] = uigetfile(fileTypeStr,failStr);
        file = [path,name];
        data = imread(file);
        filepath = file;
    end