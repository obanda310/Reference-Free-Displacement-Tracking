    function pbPillarChoice(src,~)
        disp(get(src.Parent,'selectionType'))
        if strcmp( get(src.Parent,'selectionType') , 'alt')
        src.Parent.UserData(2,1) = str2double(src.String);
        else
        src.Parent.UserData(1,1) = str2double(src.String);
        end
    end