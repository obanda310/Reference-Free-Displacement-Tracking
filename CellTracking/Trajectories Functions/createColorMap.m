function [colorMap1,colorMap2,divisionsNumber,divisionsSize,map,scheme] = createColorMap(book1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.3 Creating a Color Map for Quiver Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

promptColorMap = 'How many colors on color map? Enter an integer and press enter: ';
divisionsNumber = input(promptColorMap); %color map divisions
[map,scheme] = brewermap_view(divisionsNumber)
divisionsNumber = size(map,1);
divisionsSize = (max(max(book1(19,:,:))))/divisionsNumber;
numIndices = size(book1,1);
totalNumFrames = size(book1,2);
numTraj = size(book1,3);
colorMap1 = zeros(numIndices,totalNumFrames,numTraj,divisionsNumber); %for creating colormaps of deformations in current frame
colorMap2 = zeros(numIndices,totalNumFrames,numTraj,divisionsNumber); %for creating colormaps of largest deformation all in one frame
maskArray = zeros(numIndices,totalNumFrames,numTraj,1);
for i = 1:divisionsNumber
     colorMap1(:,:,:,i) = book1(:,:,:);     
     maskArray(1,:,:,1) = book1(19,:,:) < (divisionsSize*i) & book1(19,:,:) > (divisionsSize*(i-1));
     for j = 1:numIndices         
     colorMap1(j,:,:,i) = colorMap1(j,:,:,i).*maskArray(1,:,:,1);
     end
end
disp('done cMDBook1')
for i = 1:divisionsNumber
     colorMap2(:,:,:,i) = book1(:,:,:);
     maskArray(1,:,:,1) = book1(29,:,:) < (divisionsSize*i) & book1(29,:,:) > (divisionsSize*(i-1));
     for j = 1:numIndices         
     colorMap2(j,:,:,i) = colorMap2(j,:,:,i).*maskArray(1,:,:,1);
     end
end
disp('done cMDBook2')
end