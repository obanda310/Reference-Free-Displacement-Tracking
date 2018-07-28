function [vq2pos,vq2neg,rDisp2PosTotal, rDisp2PosMax, rDisp2NegTotal, rDisp2NegMax]  = vqStats(vq2,planesGroups)
vq2pos = vq2(:,:,:).*(vq2(:,:,:)>0);
vq2pos(isnan(vq2pos)) = 0;
vq2neg = vq2(:,:,:).*(vq2(:,:,:)<0);
vq2neg(isnan(vq2neg)) = 0;
for i = 1:size(planesGroups,1)
rDisp2PosTotal(i,1) = sum(sum(vq2pos(:,:,i)));
rDisp2PosMax(i,1) = max(max(vq2pos(:,:,i)));
rDisp2NegTotal(i,1) = sum(sum(vq2neg(:,:,i)));
rDisp2NegMax(i,1) = max(max(abs(vq2neg(:,:,i)))) * -1;
end
