function [data,dataKey] = InputSelector(auto)
    dataKey(1,1) = 2; %xCol
    dataKey(2,1) = 3; %yCol
    dataKey(3,1) = 4; %fCol
    dataKey(4,1) = 5; %tCol
    dataKey(5,1) = 7; %intCol
    dataKey(6,1) = 6; %totalCol
    if auto == 1
        w = 'Yes';
    else
    w = questdlg('Use a scalefactor of 1?',...
    'Scaling from tracking','Yes','No','Yes');
    waitfor(w);
    end
    
    
    if strcmp(w,'Yes') == 1
    dataKey(7,1) = 1;
    else
    prompt = 'What was the scale factor print out of tracking.m? Check the command window. Enter a decimal and press enter: ';
    dataKey(7,1) = input(prompt); %pixelScale
    end
    
    dataKey(8,1) = 0; %startVar
    dataKey(9,1) = .1625 * dataKey(7,1);
    dataKey(10,1) = 0.4;
    
    data = loadSpecific('trajectoriesInput.txt','*.txt','Select .txt File From Particle Tracker Output');

end