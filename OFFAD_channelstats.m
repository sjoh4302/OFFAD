channels=string(fields(OFFDATA));
for i = 1:length(fields(OFFDATA))
    for j = 1:length(fields(OFFDATA))
        coherenceMat(i,j)=length(intersect(OFFDATA.(channels(i)).AllOFFtimes,...
            OFFDATA.(channels(j)).AllOFFtimes))/length(OFFDATA.(channels(i)).AllOFFtimes);
    end
end
coherenceMat(coherenceMat==1)=NaN;
cMap = interp1([0;1],[0 1 0; 1 0 0],linspace(0,1,256));
figure
c=heatmap(coherenceMat,'Colormap',flipud(cMap));

figure
scatter(OFFDATA.(channels(1)).AllOFFtimes,...
     repmat(1,1,length(OFFDATA.(channels(1)).AllOFFtimes)))
 hold on
 scatter(OFFDATA.(channels(2)).AllOFFtimes,...
     repmat(2,1,length(OFFDATA.(channels(2)).AllOFFtimes)))
 scatter(OFFDATA.(channels(3)).AllOFFtimes,...
     repmat(3,1,length(OFFDATA.(channels(3)).AllOFFtimes)))
 scatter(OFFDATA.(channels(4)).AllOFFtimes,...
     repmat(4,1,length(OFFDATA.(channels(4)).AllOFFtimes)))
 
binEdge=[0:2:200 210:10:250 300 350]; 
figure
subplot(2,2,1)
histogram(((diff(OFFDATA.(channels(1)).nr,1,2)+1/498.2462))*1000,binEdge)
ylim([0,50])
subplot(2,2,2)
histogram(((diff(OFFDATA.(channels(2)).nr,1,2)+1/498.2462))*1000,binEdge)
ylim([0,50])
subplot(2,2,3)
histogram(((diff(OFFDATA.(channels(3)).nr,1,2)+1/498.2462))*1000,binEdge)
ylim([0,50])
subplot(2,2,4)
histogram(((diff(OFFDATA.(channels(4)).nr,1,2)+1/498.2462))*1000,binEdge)
ylim([0,50])

[median(((diff(OFFDATA.(channels(1)).nr,1,2)+1/498.2462))*1000),...
    median(((diff(OFFDATA.(channels(2)).nr,1,2)+1/498.2462))*1000)...
    median(((diff(OFFDATA.(channels(3)).nr,1,2)+1/498.2462))*1000)...
    median(((diff(OFFDATA.(channels(4)).nr,1,2)+1/498.2462))*1000)]

[mean(((diff(OFFDATA.(channels(1)).nr,1,2)+1/498.2462))*1000),...
    mean(((diff(OFFDATA.(channels(2)).nr,1,2)+1/498.2462))*1000)...
    mean(((diff(OFFDATA.(channels(3)).nr,1,2)+1/498.2462))*1000)...
    mean(((diff(OFFDATA.(channels(4)).nr,1,2)+1/498.2462))*1000)]
