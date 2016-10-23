% Required input for trajectories.m has columns:
% Count | Trajectory | Frame | x | y | z | m0 | m1 | m2 | m3 | m4 | NPscore
% This function creates columns 1 through 4

function final = createExcelForTrajectories(pillarBook)
    final = zeros(1,4);
    for a = 1:size(pillarBook,3)
        trajectory = a*ones(size(pillarBook,1),1);
        xy = pillarBook(:,1:2,a);
        frame = pillarBook(:,3,a);
        aMat = cat(2,trajectory,frame,xy);
        final = cat(1,final,aMat);
    end

    count = (1:length(final))';
    final = cat(2,count,final);
    xlswrite('trajectoriesInput.xlsx',final)
end