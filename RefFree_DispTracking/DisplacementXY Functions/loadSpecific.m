function [data,filepath] = loadSpecific(nameStr,fileTypeStr,failStr)
%Currently handles only .tif and .txt files

    strL = length(nameStr);
    files = dir(fileTypeStr); %Check Directory for default filenames
    for k = 1:length(files)
        current=files(k).name;
        if length(current)>=strL
            check(k)=strcmp(current(end-(strL-1):end),nameStr);
        end
    end
       
    loc=find(check);
    if size(loc,1)==1 && strcmp(fileTypeStr,'*.txt')==1
        data = dlmread(files(loc(1)).name);
        filepath = files(loc(1)).name;
    elseif    size(loc,1)==1 && strcmp(fileTypeStr,'*.tif')==1
        data = imread(files(loc(1)).name);  
        filepath = files(loc(1)).name;
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