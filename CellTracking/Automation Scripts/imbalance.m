%code for calculating force imbalances
function totals = imbalance(directory)
cd(directory)
load('3Ddata.mat')
shear.ltLastdX(isnan(shear.ltLastdX))= 0;
shear.ltLastdY(isnan(shear.ltLastdY))= 0;
totals(1,1) = abs(sum(shear.ltLastdX .* shear.coCheck2)); %imbalance
totals(1,3) = abs(sum(shear.ltLastdY .* shear.coCheck2)); %imbalance
totals(1,2) = sum(abs(shear.ltLastdX .* shear.coCheck2)); %total
totals(1,4) = sum(abs(shear.ltLastdY .* shear.coCheck2)); %total

zTarget = 7;
for j = 1:size(planesGroups,1)
    planesLoc3(j) = mean(planesLoc2(1,planesGroups(j,1:nnz(planesGroups(j,:)))));
end
zPlane = find(abs(planesLoc3-zTarget) == min(abs(planesLoc3-zTarget)),1,'first');
for i = 1:nnz(planesGroups(zPlane,:))
    if i == 1
        Zdisp = m3.dispFilt(1:nnz(m3.dispFilt(:,3)),3);
    else
        Zdisp = cat(1,Zdisp,m3.dispFilt(1:nnz(m3.dispFilt(:,3)),3));
    end
end
totals(1,5) = abs(sum(Zdisp)); %imbalance
totals(1,6) = sum(abs(Zdisp)); %total
totals(1,7) = ((totals(1,1)+totals(1,3))/(totals(1,2)+totals(1,4)))*100;
totals(1,8) = (totals(1,5)/totals(1,6))*100;
totals(1,9) = sum(sum(image.Area ==0))*raw.dataKey(9,1)^2;
totals(1,10) = (totals(1,2)+totals(1,4))/totals(1,6);

