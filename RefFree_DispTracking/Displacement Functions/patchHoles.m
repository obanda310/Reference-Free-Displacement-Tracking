function [r,rows,VPlanes,plane,trefAll] = patchHoles(r,rows,VPlanes,plane)
%% Attempt to repair missing data using existing incomplete entries v3
numPasses = 4;
maxDistCO = 3;
repeat = 1;
distCO = 1;
trefAll = [];
pass = 0;
progressbar('Limit','Pass')
while repeat == 1      
    progressbar(distCO/maxDistCO,pass/numPasses)
    repeatChk = 0;
    idx = find(r.vplane==0); %Grab all objects which have not been assigned to a vplane
    tplane = r.planeGroup(idx);
    idx3 = zeros(size(idx));
    cand = zeros(size(idx));
    distcheckmain = zeros(size(idx));
    alignmentbot = zeros(size(idx));
    distcheckbot = zeros(size(idx));
    alignmentleft = zeros(size(idx));
    distcheckleft = zeros(size(idx));
    alignmentright = zeros(size(idx));
    distcheckright = zeros(size(idx));
    tpmax = max(r.planeGroup);
    defCO = 0.2; %minimum deformation needed to consider alignment
    defCO2 = 0.3; %minimum deformation needed to consider alignment
    alignCO = 0.5; %minimum alignment to be considered 'aligned' (value should fall between 0 and 1, where 0 is orthogonal and 1 is parallel)
    
    for i = 1:nnz(idx) %for each unassigned object in 'idx'
        [idx2,dist2] = rangesearch([r.X,r.Y],[r.X(idx(i)),r.Y(idx(i))],4); %find nearby objects
        idx2p = idx2{1}(r.planeGroup(idx2{1}(:))==tplane(i)+1 & r.colUp(idx2{1}(:))==0 & r.vplane(idx2{1}(:))>0); %determine which of those objects are in columns that are missing data
        dist2p = dist2{1}(r.planeGroup(idx2{1}(:))==tplane(i)+1 & r.colUp(idx2{1}(:))==0 & r.vplane(idx2{1}(:))>0); %also store distance to those objects
        if nnz(idx2p)>0
            for j = 1:nnz(idx2p) % for each potential match
                %Establish reference location if current member of 'idx' is
                %matched with 'idx2p(j)'
                tref(1,1) = r.rX(idx2p(j)) + VPlanes.vX(r.vplane(idx2p(j)),tpmax-r.planeGroup(idx2p(j))+1);
                tref(1,2) = r.rY(idx2p(j)) + VPlanes.vY(r.vplane(idx2p(j)),tpmax-r.planeGroup(idx2p(j))+1);
                tref(1,3) = r.rZ(idx2p(j)) + VPlanes.vZ(r.vplane(idx2p(j)),tpmax-r.planeGroup(idx2p(j))+1);
                
                trefAll = cat(1,trefAll,tref); %store all reference locations for posterity
                
                tvec1 = [r.X(idx(i))-tref(1,1),r.Y(idx(i))-tref(1,2)]./norm([r.X(idx(i))-tref(1,1),r.Y(idx(i))-tref(1,2)]); %calculate unit displacement vector based on temporary ref 'tref'
                distcheckmain(i,j) = (norm([r.X(idx(i))-tref(1,1),r.Y(idx(i))-tref(1,2)])); % also store the magnitude
                
                %try to find bottom alignment               
                tvec2 = [r.dX(idx2p(j)),r.dY(idx2p(j))]./norm([r.dX(idx2p(j)),r.dY(idx2p(j))]); %calculate unit displacement vector for the current potential match (which will be below tvec1 if matched)
                alignmentbot(i,j) = dot((tvec1),(tvec2)./norm(tvec2))*(norm([r.dX(idx2p(j)),r.dY(idx2p(j))])>defCO); %check alignement with that vector and store if magnitude of tvec2 is greater than threshold (*alignment will be ignored if below threshold to avoid bad alignment due to noise)
                distcheckbot(i,j) = norm([r.dX(idx2p(j)),r.dY(idx2p(j))]); %store magnitude of tvec2
                
                
                %try to find left alignment
                if r.rowLeft(idx2p(j))>0 %if object to left exists
                    if r.colUp(r.rowLeft(idx2p(j)))>0 %and if object above object to left exists
                        %perform similar alignment test
                        tvec2 = [r.dX(r.colUp(r.rowLeft(idx2p(j)))),r.dY(r.colUp(r.rowLeft(idx2p(j))))]./norm([r.dX(r.colUp(r.rowLeft(idx2p(j)))),r.dY(r.colUp(r.rowLeft(idx2p(j))))]);
                        alignmentleft(i,j) = dot((tvec1),(tvec2)./norm(tvec2))*(norm([r.dX(r.colUp(r.rowLeft(idx2p(j)))),r.dY(r.colUp(r.rowLeft(idx2p(j))))])>defCO2);
                        distcheckleft(i,j) = (norm([r.dX(r.colUp(r.rowLeft(idx2p(j)))),r.dY(r.colUp(r.rowLeft(idx2p(j))))]));
                    else
                        alignmentleft(i,j) = 0;
                        distcheckleft(i,j) = 0;
                    end
                else
                    alignmentleft(i,j) = 0;
                    distcheckleft(i,j) = 0;
                end
                
                %try to find right alignment
                if r.rowRight(idx2p(j))>0 %if object to the right exists
                    if r.colUp(r.rowRight(idx2p(j)))>0 %and if object above object to the right exists
                        %perform similar alignment test
                        tvec2 = [r.dX(r.colUp(r.rowRight(idx2p(j)))),r.dY(r.colUp(r.rowRight(idx2p(j))))]./norm([r.dX(r.colUp(r.rowRight(idx2p(j)))),r.dY(r.colUp(r.rowRight(idx2p(j))))]);
                        alignmentright(i,j) = dot((tvec1),(tvec2)./norm(tvec2))*(norm([r.dX(r.colUp(r.rowRight(idx2p(j)))),r.dY(r.colUp(r.rowRight(idx2p(j))))])>defCO2);
                        distcheckright(i,j) = norm([r.dX(r.colUp(r.rowRight(idx2p(j)))),r.dY(r.colUp(r.rowRight(idx2p(j))))]);
                    else
                        alignmentright(i,j) = 0;
                        distcheckright(i,j) = 0;
                    end
                else
                    alignmentright(i,j) = 0;
                    distcheckright(i,j) = 0;
                end
            end
            
            %Convert low alignment entries to zero
            
            alignmentbot(alignmentbot>0 & alignmentbot<alignCO) = -1;
            alignmentleft(alignmentleft>0 & alignmentleft<alignCO) = -1;
            alignmentright(alignmentright>0 & alignmentright<alignCO) = -1;
            
            
            %Create 2 types of criteria for matches:
            
            %Type 1 - Alignment data exists. Three checks are performed:
            %1.) No alignment data is negative (poor alignement)
            %2.) At least 1 alignment data point is nonzero and postive
            %3.) The resulting length of 'tvec1' is greater than the
            %threshold
            tempcand = idx2p((alignmentbot(i,1:nnz(idx2p))>0 | alignmentleft(i,1:nnz(idx2p))>0 | alignmentright(i,1:nnz(idx2p))>0) & (alignmentbot(i,1:nnz(idx2p))>=0 & alignmentleft(i,1:nnz(idx2p))>=0 & alignmentright(i,1:nnz(idx2p))>=0) & distcheckmain(i,1:nnz(idx2p))>defCO2);
            tempcandd = dist2p((alignmentbot(i,1:nnz(idx2p))>0 | alignmentleft(i,1:nnz(idx2p))>0 | alignmentright(i,1:nnz(idx2p))>0) & (alignmentbot(i,1:nnz(idx2p))>=0 & alignmentleft(i,1:nnz(idx2p))>=0 & alignmentright(i,1:nnz(idx2p))>=0) & distcheckmain(i,1:nnz(idx2p))>defCO2);
            cand(i,1:nnz(tempcand)) = tempcand; %for viewing matches after
            %canddist(i,1:nnz(tempcandd)) = tempcandd; %for viewing matches after
            
            %Type 2 - All alignment tests failed displacement magnitude
            %thresholds (probably little to no real displacement here match
            %should be easy). Additionally, these should all occur on first
            %pass
            
            lowdispCand = ((distcheckbot(i,1:nnz(idx2p))<defCO) & (distcheckleft(i,1:nnz(idx2p))<defCO) & (distcheckright(i,1:nnz(idx2p))<defCO) & distcheckmain(i,1:nnz(idx2p))<defCO);
            ldcand =  idx2p(lowdispCand);
            ldcandd =  dist2p(lowdispCand);
            %ldcandAll(i,1:nnz(ldcand)) = ldcand; %for viewing matches after
            %ldcanddistAll(i,1:nnz(ldcandd)) = ldcandd; %for viewing matches after
            
            if nnz(ldcand)>0 %type 2 has priority
                repeatChk = 1;
                match = ldcand(ldcandd==min(ldcandd));                
                r.col(idx(i)) = r.col(match);
                r.vplane(idx(i)) = r.vplane(match);
                clear tempr
                tempr(:,1) = VPlanes.gridRow(r.planeGroup(idx(i)),:,r.vplane(idx(i)));
                tempr(isnan(tempr),:) = [];
                tempr(tempr==0,:) = [];
                r.row(idx(i)) = mode(tempr);
                r.colUp(match)= idx(i); %Need to update this so the column is not matched more than once
                
                %Update displacement data for use in future matches'
                %alignment
                r.rX(idx(i)) = r.rX(match) + VPlanes.vX(r.vplane(match),tpmax-r.planeGroup(match)+1);
                r.rY(idx(i)) = r.rY(match) + VPlanes.vY(r.vplane(match),tpmax-r.planeGroup(match)+1);
                r.rZ(idx(i)) = r.rZ(match) + VPlanes.vZ(r.vplane(match),tpmax-r.planeGroup(match)+1);
                r.dX(idx(i)) = r.X(idx(i))-r.rX(idx(i));
                r.dY(idx(i)) = r.Y(idx(i))-r.rY(idx(i));
                r.dZ(idx(i)) = r.Z(idx(i))-r.rZ(idx(i));
                r.dS(idx(i)) = sqrt(r.dX(idx(i))^2 + r.dY(idx(i))^2);
                
            elseif min(tempcandd)<distCO %then check for type 1
                repeatChk = 1;
                match = tempcand(tempcandd==min(tempcandd));
                r.col(idx(i)) = r.col(match);
                r.vplane(idx(i)) = r.vplane(match);
                clear tempr
                tempr(:,1) = VPlanes.gridRow(r.planeGroup(idx(i)),:,r.vplane(idx(i)));
                tempr(isnan(tempr),:) = [];
                tempr(tempr==0,:) = [];
                r.row(idx(i)) = mode(tempr);
                r.colUp(match)= idx(i); %Need to update this so the column is not matched more than once
                
                %Update displacement data for use in future matches'
                %alignment
                r.rX(idx(i)) = r.rX(match) + VPlanes.vX(r.vplane(match),tpmax-r.planeGroup(match)+1);
                r.rY(idx(i)) = r.rY(match) + VPlanes.vY(r.vplane(match),tpmax-r.planeGroup(match)+1);
                r.rZ(idx(i)) = r.rZ(match) + VPlanes.vZ(r.vplane(match),tpmax-r.planeGroup(match)+1);
                r.dX(idx(i)) = r.X(idx(i))-r.rX(idx(i));
                r.dY(idx(i)) = r.Y(idx(i))-r.rY(idx(i));
                r.dZ(idx(i)) = r.Z(idx(i))-r.rZ(idx(i));
                r.dS(idx(i)) = sqrt(r.dX(idx(i))^2 + r.dY(idx(i))^2);
                
            end
        end
        idx3(i,1:nnz(idx2p)) = idx2p;
        %dist3(i,1:nnz(dist2p)) = dist2p;
    end
    
    if repeatChk == 1
        repeat = 1;
        pass= pass+1;
        if pass == numPasses
            distCO = distCO+1;
            pass = 0;
        end
    elseif distCO <maxDistCO
        repeat = 1;
        distCO = distCO+1;
    else
        repeat = 0;
    end
    trefAll = unique(trefAll,'rows');
    
    % Clean and update all major class variables before repeating
    if repeatChk == 1
        [r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);
    end   
    
    if distCO>maxDistCO
        repeat = 0;
    end

end
%% Auto-Correct Cont.
%Grab the average displacement in a small radius, apply that to missing
%location, and determine if a match exists in the vicinity, if yes, match
%to that one. 
%****Currently only operates on top 2 planes, may need some
%generalization in future****

for i = 1:3
tpmax = max(r.planeGroup);
idx2p = find(r.colUp == 0 & r.planeGroup==2 & r.vplane>0);
tidx = find(r.planeGroup==1);
topPlane = find(r.planeGroup==1);
t = delaunayn(double(r.r(topPlane,1:2)));
trefAll = [];
pidx2 = [];
for j = 1:nnz(idx2p) % for each potential match
    %Establish reference location if current member of 'idx' is
    %matched with 'idx2p(j)'
    tref(1,1) = r.rX(idx2p(j)) + VPlanes.vX(r.vplane(idx2p(j)),tpmax-r.planeGroup(idx2p(j))+1);
    tref(1,2) = r.rY(idx2p(j)) + VPlanes.vY(r.vplane(idx2p(j)),tpmax-r.planeGroup(idx2p(j))+1);
    tref(1,3) = r.rZ(idx2p(j)) + VPlanes.vZ(r.vplane(idx2p(j)),tpmax-r.planeGroup(idx2p(j))+1);
    trefAll = cat(1,trefAll,tref); %store all reference locations for posterity
    [pidx,~] = rangesearch([r.rX(tidx) r.rY(tidx)],tref(1,1:2),5); %range was 3 before instead of 5
    pidx2 = cat(1,pidx2,tidx(pidx{1}(:)));
    meanV = mean([r.dX(tidx(pidx{1}(:))),r.dY(tidx(pidx{1}(:)))],1);
    meanV2{j} = mean([r.dX(tidx(pidx{1}(:))),r.dY(tidx(pidx{1}(:)))],1);
    meanV2{j}(1,3) = 0;
    [pidx,pdist] = dsearchn(double(r.r(topPlane,1:2)),t,tref(1,1:2)+meanV);
    if pdist < 1 && (r.row(topPlane(pidx))>0||r.col(topPlane(pidx))>0) % if its a good match, but feature is already linked
        if r.colDown(topPlane(pidx))>0
        r.colUp(r.colDown(topPlane(pidx))) = 0; %unlink previous pillar
        end
        r.row(topPlane(pidx)) = mode(VPlanes.gridRow(r.planeGroup(topPlane(pidx)),:,r.vplane(idx2p(j))));
        r.col(topPlane(pidx)) = r.col(idx2p(j));
        r.colDown(topPlane(pidx)) = idx2p(j);        
    elseif pdist < 1.5 && (r.row(topPlane(pidx))== 0 && r.col(topPlane(pidx))==0)  %otherwise, if it is an okay match ****there may be a bug in using (r.row(topPlane(pidx))== 0)**** the zero should probably be any row not associated with the current incomplete column's vertical plane            
        r.row(topPlane(pidx)) = mode(VPlanes.gridRow(r.planeGroup(topPlane(pidx)),:,r.vplane(idx2p(j))));
        r.col(topPlane(pidx)) = r.col(idx2p(j));               
        r.colDown(topPlane(pidx)) = idx2p(j);
    end         
end
% Debugging
% ViewVPlanes(VPlanes,r)
% meanV3 = cat(1,meanV2{:});
% scatter3(r.X((pidx2)),r.Y((pidx2)),r.Z((pidx2)))
% scatter3(trefAll(:,1),trefAll(:,2),trefAll(:,3))
% quiver3(trefAll(:,1),trefAll(:,2),trefAll(:,3),meanV3(:,1),meanV3(:,2),meanV3(:,3),0)
[r,VPlanes,rows,plane] = RefreshVars(r,rows,plane);
end

trefAll = unique(trefAll,'rows');