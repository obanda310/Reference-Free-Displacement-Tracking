function AR = cellAspectRatio(directory)
if nargin ==1
    cd(directory);
end
close all

%%
load('3Ddata.mat')
props = regionprops(image.Area==0,'MajorAxisLength','MinorAxisLength');
Ma = cat(1,(props.MajorAxisLength));
MaMax = max(Ma);
Mi = cat(1,(props.MinorAxisLength));
MiMax = max(Mi);

AR(1,1) = MaMax/MiMax;
AR(1,2) = MaMax*raw.dataKey(9,1);
AR(1,3) = MiMax*raw.dataKey(9,1);