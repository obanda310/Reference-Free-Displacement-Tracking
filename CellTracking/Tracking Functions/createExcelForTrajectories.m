function final = createExcelForTrajectories(pillarBook)
    final = zeros(1,7);
    for a = 1:size(pillarBook,3)
        first = find(pillarBook(:,1,a),1,'first'); %these two lines will get rid of zeros in the output excel sheet
        last = find(pillarBook(:,1,a),1,'last');   %
        xy = pillarBook(first:last,1:2,a);
        frame = pillarBook(first:last,3,a);
        trajectory = pillarBook(first:last,4,a);
        intensity = pillarBook(first:last,5,a);
        intensityGauss = pillarBook(first:last,6,a);
        intensityMean = pillarBook(first:last,7,a);
        aMat = cat(2,xy,frame,trajectory,intensity,intensityGauss,intensityMean);
        final = cat(1,final,aMat);
    end

    count = (1:length(final))';
    final = cat(2,count,final);
    
    excelFileName = 'trajectoriesInput.txt';
    % If a file already exists with the name of the Excel spreadsheet 
    % specified in excelFileName, that file is deleted. This is necessary 
    % because the file will not simply be overwritten. If the array being
    % saved currently is smaller in either dimension than the array saved
    % previously, the elements of the previous array that are outside the
    % bounds of the current array will remain, and not be deleted or
    % overwritten.
    if exist(excelFileName,'file')
        delete(excelFileName)
    end
    dlmwrite(excelFileName,final);
end