function MoveMatFiles(directory)
if nargin == 1
    cd(directory)
end

files = dir('*.mat');
mkdir('Matlab Data Files')
for i = 1:size(files,1)
    movefile(files(i).name,'Matlab Data Files')
end

try
movefile('Shear Mat Files\DataRaw.mat','Matlab Data Files\DataRaw.mat')
catch
end