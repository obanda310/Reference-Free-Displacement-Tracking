function [colorMap1,colorMap2,divisionsNumber,divisionsSize,map,scheme] = createColorMap(book1,book2,outputs)

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

divisionsSize = 1;
numIndices = size(book1,1);
totalNumFrames = size(book1,2);
numTraj = size(book1,3);
%colorMap1 = zeros(numIndices,totalNumFrames,numTraj,divisionsNumber); %for creating colormaps of deformations in current frame
%colorMap2 = zeros(numTraj,size(book2,2),divisionsNumber); %for creating colormaps of largest deformation all in one frame
%maskArray1 = zeros(numIndices,totalNumFrames,numTraj,1);
%maskArray2 = zeros(numTraj,1,1);
for i = 1:divisionsNumber
    for j = 1:size(book1,2)
    clear maskArray
     if i == divisionsNumber
     maskArray1 = book1(15,j,:) > (divisionsSize*(i-1));
     else
     maskArray1 = book1(15,j,:) <= (divisionsSize*i) & book1(15,j,:) > (divisionsSize*(i-1));
     end         
     colorMap1(1:size(find(squeeze(maskArray1)),1),i,j) = find(squeeze(maskArray1));
    end
end
for i = 1:divisionsNumber
    clear maskArray
     if i == divisionsNumber
     maskArray1 = book2(:,20) > (divisionsSize*(i-1));
     else
     maskArray1 = book2(:,20) <= (divisionsSize*i) & book2(:,20) > (divisionsSize*(i-1));
     end         
     colorMap2(1:size(find(maskArray1(:,:)),1),i) = find(maskArray1(:,:));
end
disp('done cMDBook1')
% for i = 1:divisionsNumber
%      colorMap2(:,:,i) = book2(:,:,:);
%      if i == divisionsNumber
%           maskArray2(:,1,1) = book2(:,9) > (divisionsSize*(i-1));
%      else
%      maskArray2(:,1,1) = book2(:,9) <= (divisionsSize*i) & book2(:,9) > (divisionsSize*(i-1));
%      end
%      for j = 1:numIndices         
%      colorMap2(:,j,i) = colorMap2(:,j,i).*maskArray2(:,1,1);
%      end
% end
% disp('done cMDBook2')

end