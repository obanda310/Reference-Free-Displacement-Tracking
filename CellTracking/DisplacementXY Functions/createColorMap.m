function [colorMap1,colorMap2,divisionsNumber,divisionsSize,map,scheme] = createColorMap(shear,raw,outputs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.2 Creating a Color Map for Quiver Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set up color map
if ismember(1,outputs) == 1
    divisionsNumber = 30; % default value
    map = brewermap(divisionsNumber,'*Spectral');
    scheme = 'Spectral';
else
promptColorMap = 'How many colors on color map? Enter an integer and press enter: ';
divisionsNumber = input(promptColorMap); %color map divisions
[map,scheme] = brewermap_view(divisionsNumber);
divisionsNumber = size(map,1);
end
divisionsSize = 1*raw.dataKey(9,1);

for i = 1:divisionsNumber
    for j = 1:shear.numFrames
    clear maskArray
     if i == divisionsNumber
     maskArray1 = shear.ltdXY(j,:) >= (divisionsSize*(i-1));
     else
     maskArray1 = shear.ltdXY(j,:) < (divisionsSize*i) & shear.ltdXY(j,:) >= (divisionsSize*(i-1));
     end 
     
     colorMap1(1:size(find(squeeze(maskArray1)'),1),i,j) = find(squeeze(maskArray1))';
    end
end

for i = 1:divisionsNumber
    clear maskArray
     if i == divisionsNumber
     maskArray1 = shear.ltLastdXY(:) > (divisionsSize*(i-1));
     else
     maskArray1 = shear.ltLastdXY(:) <= (divisionsSize*i) & shear.ltLastdXY(:) > (divisionsSize*(i-1));
     end         
     colorMap2(1:size(find(maskArray1(:,:)),1),i) = find(maskArray1(:,:))';
end

disp('done cMDBook1')

end