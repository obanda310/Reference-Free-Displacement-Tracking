
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 3.) Identifying Neighborhood Trajectories for Each Trajectory
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% disp('3.1 Identifying Pillar Neighbor Groups')
% 
% clear tempDistances
% book3 = zeros(numTraj,50);
% cm3 = zeros(numTraj,50,totalNumFrames);
% cm4 = zeros(numTraj,totalNumFrames);
% interpXq = linspace(1,size(book1,2),size(book1,2)*3);
% maxDistance = (max(max(book1(5,:,:))));
% for i = 1:numTraj
%     tempDistances(1:numTraj,1) = ((book2(:,1)-book2(i,1)).^2.+(book2(:,2)-book2(i,2)).^2).^0.5;
%     tempDistances(1:numTraj,2) = linspace(1,numTraj,numTraj);
%     tempDistances = sortrows(tempDistances);
%     book3(i,1:50) = tempDistances(1:50,2); %record the closest 50 pillars
%     if book2(i,9) < (maxDistance) %filter by deformation limits here and in following if statement
%         cm4(i,1:totalNumFrames) = book1(6,1:totalNumFrames,i);
%     else
%         cm4(i,1:totalNumFrames) = NaN;
%     end
%     
%     for j = 1:50
%         if book2(book3(i,j),9) < (maxDistance) % *following if statement*
%             cm3(i,j,1:totalNumFrames) = book1(6,1:totalNumFrames,book3(i,j)); %record intensity value if pillar has low deformation
%         else
%             cm3(i,j,1:totalNumFrames) = NaN;
%         end
%     end
% end
% 
% cm4(cm4==0) = NaN;
% totalAverage = zeros(1,totalNumFrames);
% for i = 1:totalNumFrames
%     totalAverage(1,i) = mean(cm4(:,i),'omitnan');
% end
% totalAverageInterp(1,:) = interp1(linspace(1,size(book1,2),size(book1,2)),conv(totalAverage(1,:),gausswin(6),'same'),interpXq,'spline');
% 
