% Required input for trajectories.m has columns:
% Count | Trajectory | Frame | x | y | z | m0 | m1 | m2 | m3 | m4 | NPscore
% This function creates columns 1 through 4

function final = createExcelForTrajectories(pillarBook)
    final = zeros(1,5);
    for a = 1:size(pillarBook,3)
        first = find(pillarBook(:,1,a),1,'first'); %these two lines will get rid of zeros in the output excel sheet
        last = find(pillarBook(:,1,a),1,'last');   %
        xy = pillarBook(first:last,1:2,a);
        frame = pillarBook(first:last,3,a);
        trajectory = pillarBook(first:last,4,a);
        intensity = pillarBook(first:last,5,a);
        aMat = cat(2,xy,frame,trajectory,intensity);
        final = cat(1,final,aMat);
    end

    count = (1:length(final))';
    final = cat(2,count,final);
    xlswrite('trajectoriesInput.xlsx',final);
end