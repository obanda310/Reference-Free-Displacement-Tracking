function clearDirectory(directory)
if nargin ==1
    cd(directory)
end

% try
%     rmdir('3D Plots','S')
% catch
% end
% 
% try
%    rmdir('HeatMaps','S')
% catch
% end
% 
% try
%    rmdir('Histograms','S')
% catch
% end

try
    rmdir('ColorBar','S')
catch
end


try
    rmdir('Profile Data','S')
catch
end


try
    rmdir('Shear Mat Files','S')
catch
end