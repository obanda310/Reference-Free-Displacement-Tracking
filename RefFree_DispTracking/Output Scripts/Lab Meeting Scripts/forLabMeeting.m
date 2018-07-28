figure
scatter3(fullData2(:,4),fullData2(:,5),fullData2(:,6))

figure
scatter3(fDzp(:,1),fDzp(:,2),fDzp(:,3))


%%
figure
quiver3(fDz(:,1),fDz(:,2),fDz(:,3),fDzp(:,1)-fDz(:,1),fDzp(:,2)-fDz(:,2),fDzp(:,3)-fDz(:,3),0)

%%
figure
quiver3(m3.refSC(:,1),m3.refSC(:,2),m3.refSC(:,3),m3.dispFilt(:,1),m3.dispFilt(:,2),m3.dispFilt(:,3),0)
hold on
quiver3(em3.refSC(:,1),em3.refSC(:,2),em3.refSC(:,3),em3.dispFilt(:,1),em3.dispFilt(:,2),em3.dispFilt(:,3),0)
plot(Surface2)

%%
figure
quiver3(m3.refSC(:,1),m3.refSC(:,2),m3.refSC(:,3),m3.dispFilt(:,1),m3.dispFilt(:,2),m3.dispFilt(:,3),0)
hold on
quiver3(em3.refSC(:,1),em3.refSC(:,2),em3.refSC(:,3),em3.dispFilt(:,1),em3.dispFilt(:,2),em3.dispFilt(:,3),0)

quiver3(dZerosZ(:,1),dZerosZ(:,2),dZerosZ(:,3),dZerosZ(:,4),dZerosZ(:,5),dZerosZ(:,6),0)
