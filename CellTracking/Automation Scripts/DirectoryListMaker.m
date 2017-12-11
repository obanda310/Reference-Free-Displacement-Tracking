% This script should allow for creating a list of directories for automation of tracking on multiple image stacks
%% Run this section to start a new file
clear dirList
dirList =  {0};
[ListName,ListPath] = uiputfile('*.*','Open a list of Directories','List of Directories.mat');
fullListName = strcat(ListPath,ListName);
save(fullListName,'dirList')
%% Run this section to select a previously selected file
[ListName,ListPath] = uigetfile;
fullListName = strcat(ListPath,ListName);
load(strcat(ListPath,ListName),'dirList')
%% Run this section to add the current directory
AddPath = cd;
if size(dirList,1) == 0
    dirList{1,1} = AddPath
elseif dirList{1,1} == 0
    dirList{1,1} = AddPath
else
    dirList = cat(1,dirList,AddPath)
end
save(fullListName,'dirList')
%% Run this section to remove the last entry
dirList(end,:) = []


%% Run This Section to remove a specific entry
remove = [1 2 4 5 15 16]; %choose which one to remove here
dirList(remove,:) = []
%%
[ListName,ListPath] = uiputfile('*.*','Open a list of Directories','List of Directories.mat');
fullListName = strcat(ListPath,ListName);
save(fullListName,'dirList')
%% Run this section to browse and add directory
[AddName,AddPath] = uigetdir;


%% Concatenate previous directories to a new one
clear dirList
dirList =  {0};
parts = strsplit(pwd, '\');
prefix = parts{end};
[ListName,ListPath] = uiputfile('*.*','Choose Save Location',strcat(prefix,' List of Directories.mat'));

fullListName = strcat(ListPath,ListName);
save(fullListName,'dirList')
num2cat = input('How many directory lists are being concatenated?');
for i = 1:num2cat
    if i == 1
        [ListName,ListPath] = uigetfile('*.mat','Select First List');
        newdir = load(strcat(ListPath,ListName),'dirList');
        dirList = newdir.dirList;
    else
        [ListName,ListPath] = uigetfile('*.mat','Select Another List');
        adddir = load(strcat(ListPath,ListName),'dirList');
        addDir = adddir.dirList;
        dirList = cat(1,dirList,addDir);
    end
end
save(fullListName,'dirList')

