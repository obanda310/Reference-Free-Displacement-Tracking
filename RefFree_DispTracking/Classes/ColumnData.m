classdef ColumnData
    properties
        m
        
        
    end
    methods
        function obj = ColumnData(r)
            obj.m = zeros(max(r.col),max(r.plane));
            for i = 1:size(r.col,1) 

                if r.col(i,1)>0
                obj.m(r.col(i,1),r.plane(i,1)) = i;
                end
            end
        end        
    end
end