function [r] = match2D(r,raw,shear,rowsNDCU)

progressbar('Matching 3D Detections to 2D Pillars')

%r.col = zeros(1,r.l);
tempZ = r.Z;
tempcol = zeros(r.l,2);
for i = 1:r.l
    progressbar(i/r.l)
    tl = floor(tempZ(i)/raw.dataKey(10,1))-1;
    if tl <=0
        tl=1;
    end
    
    tu = ceil(tempZ(i)/raw.dataKey(10,1))+1;
    if tu <= 0
        tu = 2;
    end
    if tu > size(shear.rawX,1)
        tu = size(shear.rawX,1);
    end
    if tl > size(shear.rawX,1)
        tl = size(shear.rawX,1);
    end

    
    differences = min(squeeze(sqrt((shear.rawX(tl:tu,:)-r.X(i)).^2+(shear.rawY(tl:tu,:)-r.Y(i)).^2+(shear.rawZ(tl:tu,:)-r.Z(i)).^2)));
    if min(differences)<.5
        tempcol(i,1) = find(differences==min(differences));
        differences2 = squeeze(sqrt((shear.rawX(tl:tu,tempcol(i,1))-r.X(i)).^2+(shear.rawY(tl:tu,tempcol(i,1))-r.Y(i)).^2+(shear.rawZ(tl:tu,tempcol(i,1))-r.Z(i)).^2));
        tempcol(i,2) = tl + find(differences2==min(differences2)) - 1;
        
    end
    
end

r.col = tempcol;
disp(['done Matching 3D detections to 2D-based pillars at ' num2str(toc) ' seconds'])

%% Row Shift Correction
% figure
% hold on
% for i=1:size(rows,1)
% scatter3(r.X(rowsNDCU(i,(r.col(rowsNDCU(i,rowsNDCU(i,:)>0))>0),:)),r.Y(rowsNDCU(i,(r.col(rowsNDCU(i,rowsNDCU(i,:)>0))>0),:)),r.Z(rowsNDCU(i,(r.col(rowsNDCU(i,rowsNDCU(i,:)>0))>0),:)))
% end
%
for i = 1:size(rowsNDCU,1)
    %rowsSCtest(i,1:length(r.X(rowsNDCU(i,(r.col(rowsNDCU(i,rowsNDCU(i,:)>0))>0),:)))) = r.X(rowsNDCU(i,(r.col(rowsNDCU(i,rowsNDCU(i,:)>0))>0),:));
    rowsSC(i,1) = mean(r.X(rowsNDCU(i,(r.col(rowsNDCU(i,rowsNDCU(i,:)>0),1)>0),:))-shear.rawX1(r.col(rowsNDCU(i,(r.col(rowsNDCU(i,rowsNDCU(i,:)>0),1)>0),:),1))');
    %rowsSC(i,3) = length(r.X(r.col(rowsNDCU(rowsNDCU(i,:)>0))>0)-shear.rawX1(r.col(r.col(rowsNDCU(rowsNDCU(i,:)>0))>0))');
    rowsSC(i,2) = mean(r.Y(rowsNDCU(i,(r.col(rowsNDCU(i,rowsNDCU(i,:)>0),1)>0),:))-shear.rawY1(r.col(rowsNDCU(i,(r.col(rowsNDCU(i,rowsNDCU(i,:)>0),1)>0),:),1))');
end

%store shift correction in r
for i = 1:r.l
    r.XSC(i) = rowsSC(r.row(i),1);
    r.YSC(i) = rowsSC(r.row(i),2);
end