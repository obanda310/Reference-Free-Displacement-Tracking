function [obj,r] = buildRowsTest(obj,r,plane,preview)
%
%Note: Currently if there are seemingly duplicate features (two detections per ellipsoid),
%the code will break a row at that location. Stricter feature detection can
%solve this.
r.row = zeros(size(r.X));
obj.m = [];


rowN = 0;
%
for pass = 1:2 %first pass: assign all non-deformed objects to rows; second pass: try to assign everthing else
    for j = 1:r.l % run through every object
        if pass == 1 && ismember(j,rlist) && r.row(j) == 0 %check if current object has been assigned
            %set limits for row filters (microns)
            limit1 = 1; %distance from line (rowV passing through j)            
            rlist = r.ND; %subtract from this list to remove used points
            rLUCand = [];
            rowN = rowN+1;
            [jrow,jplane] = find(plane.final==j);
            jplane = unique(jplane);
            rlistP = plane.final(1:nnz(plane.final(:,jplane)),jplane);
            clear currentRow
            r.row(j) = rowN; %store row number
            rlist(rlist==j) = 0; % remove j from list
            clear rLU
            rLU = intersect(rlist,rlistP); %grab all non-assigned objects in the same plane
            rLU(rLU==0,:) = [];
            if size(rLU,1) == 0
                clear rLU
                rLU = 0;
            end
            clear differences dv22
            
            %narrow potential row members with limit 1
            for k=1:size(rLU,1)
                if rLU(k,1) >0
                    tempd = zeros(1,3);
                    tempd(1,1:2) = r.r(rLU(k,1),1:2)-r.r(j,1:2);
                    differences(k,1) = norm(cross(obj.V,tempd))/norm(obj.V);
                    
                else
                    differences(k,1)= 1000;
                end
            end
            rLUCand = rLU(differences(:,1)<limit1); %final list of potential matches
            r.row(rLUCand) = rowN;
            
            
            %         else % pass == 2, matching deformed objects now
            %         %set limits for row filters (microns)
            %         limit1 = 2; %distance from line (rowV passing through j)
            %         limit2 = 1; %distance from point (rowV+j)
            %         rlist = r.D; %subtract from this list to remove used points
            %             %rowN = r.row(j);
            %             rLUCand = [];
            %
            %             %PART 2
            %             %build rows feature by feature, assumption is that the closest
            %             %feature to j+/-rowV is the next row member
            %             if size(rLUCand,1)>0
            %                 %grab cross products for trajectory information
            %                 traj = zeros(size(rLUCand,1),3);
            %                 for cp = 1:size(rLUCand,1)
            %                     tempd(1,1:2) = r.r(rLUCand(cp,1),1:2)-r.r(j,1:2);
            %                     traj(cp,1:3) = cross(obj.V,tempd);
            %                 end
            %                 %direction1
            %                 [obj,r,rlist] = growRow(j,r,obj,rLUCand,traj,rlist,1,limit2,rowN);
            %                 %direction2
            %                 [obj,r,rlist] = growRow(j,r,obj,rLUCand,traj,rlist,-1,limit2,rowN);
            %             end
        end
    end
end


for i = 1:rowN
    obj.m(i,1:size(find(r.row(:)==i))) = find(r.row(:)==i);
end

if preview > 0
    figure
    hold on
    for i = 1:size(obj.m,1)
        scatter3(r.X(obj.m(i,1:nnz(obj.m(i,:)))),r.Y(obj.m(i,1:nnz(obj.m(i,:)))),r.Z(obj.m(i,1:nnz(obj.m(i,:)))))
    end
end

    function [rows,r,RemainList] = growRow(ridx,r,rows,CandList,traj,RemainList,direction,limit,currentRow)
        trajh = [0,0,0]; % this variable will hold traj history for current row
        dist(:,1) = sqrt(((r.X(ridx)*ones(size(CandList,1),1))+rows.V(1,1)-r.X(CandList(:,1))).^2+((r.Y(ridx)*ones(size(CandList,1),1))+rows.V(2,1)-r.Y(CandList(:,1))).^2);
        d1check = find(dist<limit);
        if size(d1check,1)>0
            d1proceed = 1;
            seed = ridx; %this is the index of the last rowsect added to a row.
            rowVS = direction; %this is a multiplier for the row vector that increases if no new row member is found
            
            while d1proceed == 1
                clear differences
                dist(:,1) = sqrt(((r.X(seed)*ones(size(CandList,1),1))+(rows.V(1,1)*rowVS)-r.X(CandList(:,1))).^2+(((r.Y(seed)*ones(size(CandList,1),1))+(rows.V(2,1)*rowVS))-r.Y(CandList(:,1))).^2);
                rLUMatch = CandList(dist<limit);
                
                if size(rLUMatch,1)>0
                    seed = CandList((dist==min(dist)));
                    trajh(end+1,1:3) = traj((dist==min(dist)),1:3);
                    r.row(seed,1) = currentRow;
                    RemainList(RemainList(:,1)==seed,:) = [];
                    rowVS = direction; %reset rowVS if match is found
                    
                else
                    rowVS =rowVS+direction; %increase magnitude of rowVS if no match is found
                    if r.X(seed)+(rows.V(1,1)*rowVS)>max(r.X(:))+rows.VL || r.Y(seed)+(rows.V(2,1)*rowVS)>max(r.Y(:))+rows.VL || r.X(seed)+(rows.V(1,1)*rowVS)<min(r.X(:))-rows.VL || r.Y(seed)+(rows.V(2,1)*rowVS)<min(r.Y(:))-rows.VL
                        d1proceed=0;
                    end
                end
                
            end
        end
    end
end