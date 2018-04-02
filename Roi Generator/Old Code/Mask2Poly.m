clear all
close all
IM = getImages();
%%
clear A B C IM2
close all
div = 10;
IMsize = size(IM,1);
IMdiv = round(IMsize/div);
A= ones(IMsize,IMsize);
D= zeros(IMsize,IMsize);
for i =1:div-1
A(IMdiv*i-1:IMdiv*i+1,:) = 0;   
%A(:,IMdiv*i) = 0;
end
IM2 = (IM.*A);
C = bwboundaries(IM2,4);
for i = 1:size(C,1)
   thisbound = C{i,1};
   for j = 1:size(thisbound,1)
   D(thisbound(j,1),thisbound(j,2)) = 1;
   end
end
B = bwlabel(IM2);
%corners = detectHarrisFeatures(IM,'FilterSize',41);%,
corners = detectMinEigenFeatures(IM,'FilterSize',3);
imshow(D,[])
figure
imshow(IM,[])
hold on;
plot(corners.selectStrongest(size(corners,1)));%

    h = 1/3*ones(3,1);
    H = h*h';
    % im be your image
    imfilt = filter2(H,D);
figure
se = strel('disk',10)
imshow(imdilate(imfilt>.57,se)+D,[])
%%
polys = mask2poly(IM2);
% %%
% colors=['b' 'g' 'r' 'c' 'm' 'y'];
% 
% for k=1:length(C)
%     perim = C{k}; %Obtains perimeter coordinates (as a 2D matrix) from the cell array
%     cidx = mod(k,length(colors))+1;% Obtains colours for the plot
%     plot(perim(:,2), perim(:,1),...
%         colors(cidx),'LineWidth',2);
%     
% end
% 
% Coordmat = cell2mat(C) %Converts the traced regions to a matrix
% 
% X = Coordmat(:,1)
% Y = Coordmat(:,2) % This gives the edge coordinates in matrix form
% 
% %% Corner Finding Section (from Jonas' answer to a previous question
% %# get corners
% cornerProbability = cornermetric(IM2);
% 
% 
% cornerIdx = find(cornerProbability==max(cornerProbability(:)));
% 
% %# Label the image. bwlabel puts 1 for the first feature, 2 for the second, etc.
% %# Since concave corners are placed just outside the feature, grow the features
% %# a little before labeling
% bw2 = imdilate(IM2,ones(3));
% labeledImage = bwlabel(bw2);
% 
% %# read the feature number associated with the corner
% cornerLabels = labeledImage(cornerIdx);
% 
% %# find all corners that are associated with feature 1
% corners_1 = cornerIdx(cornerLabels==1)
% 
% [Xcorners, Ycorners] = ind2sub(size(IM2),corners_1) % Convert subscripts
% 
% % Here comes the new part
% 
% %# count corners
% nCorners = length(Xcorners);
% %# corner2boundaryIdx will contain the index of the boundary pixel
% %# corresponding to a given corner
% corner2boundaryIdx = zeros(nCorners,1);
% 
% %# remove last entry of Coordmat, because the first and the last pixel are
% %# repeated
% Coordmat(end,:) = [];
% 
% %# loop since we need to tread convex and concave corners differently
% for i = 1:nCorners
%     %# find distance of corner to boundary
%     dist = sqrt((Coordmat(:,1)-Xcorners(i)).^2 + (Coordmat(:,2)-Ycorners(i)).^2);
%     %# find index of closest boundary pixel. Use find
%     %# instead of 2nd output of min, because we need to know how
%     %# many closest indices there are in order to distinguish
%     %# concave and convex corners. Convex corners are directly
%     %# on the boundary, thus there is one minimum distance (0).
%     %# Concave corners are just outside the boundary, thus, there
%     %# are two minimum distances (1)
%     minDistIdx = find(dist == min(dist));
%     
%     %# count how many minimum distances we have
%     switch length(minDistIdx)
%         case 1
%             %# convex corners. Everything is simple
%             corner2boundaryIdx(i) = minDistIdx;
%         case 2
%             %# concave corners. Take the index right in the middle. Note
%             %# for this to work, you need to have 4-connected boundaries,
%             %# otherwise, the concave corner is not part of the boundary
%             %# becasue the boundary jumps along the diagonal.
%             %# If you have to use 8-connected boundaries, the 'good'
%             %# difference is 1, and you need to calculate the coordinate of
%             %# the concave corner from the "corner" and the two adjacent
%             %# boundary pixels.
%             if diff(minDistIdx) == 2
%                 corner2boundaryIdx(i) = mean(minDistIdx);
%             else
%                 error('boundary starts/stops at concave corner - case not handled yet')
%             end
%         otherwise
%             error('%i minDist found for corner # %i',length(minDistIdx),i)
%     end
% end
% 
% %# All that is left to do is read the boundary pixel coordinates in the
% %# right order
% corner2boundaryIdx = sort(corner2boundaryIdx);
% orderedCorners = Coordmat(corner2boundaryIdx,:);
% 
% %# plot corner numbers for confirmation
% hold on
% for i = 1:nCorners
%     text(orderedCorners(i,2),orderedCorners(i,1),num2str(i),'Color','r')
% end