function [OFFDATA]=OFFAD_channelstats(OFFDATA,OFFDATA_var)

%Close loading window
g.Cluster = findobj('tag', 'OFFAD_CLUSTER');
set(g.Cluster,'Visible','off')

OFFDATA_var.ChannelIDs=string(fields(OFFDATA));
OFFDATA_var.ChannelNumsOnly=double(string(cellfun(@(X) regexp(X,'\d*','match'),OFFDATA_var.ChannelIDs,'UniformOutput',false)));

%Plotting
g.Channelstats = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576/2 432], ...
	'PaperUnits','points', ...
	'name', 'OFFAD (OFF_period Automated Detection)', ... 
	'numbertitle', 'off', ...
	'Position',[200 100 800 450], ...
    'Tag','OFFAD_CHANNELSTATS');

%%%%%%% Off durations plot (default)
OFFdurations=[];
OFFdurationsID=[];
for i = 1:length(OFFDATA_var.ChannelIDs)
   OFFdurations=[OFFdurations; (diff(OFFDATA.(OFFDATA_var.ChannelIDs(i)).nr,1,2)+1/498.2462)*1000];
   OFFdurationsID=[OFFdurationsID; repmat(OFFDATA_var.ChannelNumsOnly(i),length(OFFDATA.(OFFDATA_var.ChannelIDs(i)).nr),1)];
end
subplot('position',[0.07 0.12 0.65 0.8])
violinplot(OFFdurations,OFFdurationsID,...
'ShowData',false,'ViolinAlpha',0.3);
ylabel('Off period duration (ms)');

 
%%%%%%% Channel coherence plot 
for i = 1:length(fields(OFFDATA))
    for j = 1:length(fields(OFFDATA))
        coherenceMat(i,j)=length(intersect(OFFDATA.(OFFDATA_var.ChannelIDs(i)).AllOFFtimes,...
            OFFDATA.(OFFDATA_var.ChannelIDs(j)).AllOFFtimes))/length(OFFDATA.(OFFDATA_var.ChannelIDs(j)).AllOFFtimes);
    end
end
coherenceMat(coherenceMat==1)=NaN;
coherence=reshape(coherenceMat,[],1);
coherenceID=reshape(repmat(OFFDATA_var.ChannelNumsOnly',length(OFFDATA_var.ChannelNumsOnly),1),[],1);



%%%%%%%%% Off period number
OFFpnum=[];
currentFields=fields(OFFDATA.(OFFDATA_var.ChannelIDs(1)));
for i = 1:length(OFFDATA_var.ChannelIDs)
     OFFpnum=[OFFpnum; length(OFFDATA.(OFFDATA_var.ChannelIDs(i)).(string(currentFields(1))))];
end



%%%%%%%%% Off period number
OFFtime=[];
currentFields=fields(OFFDATA.(OFFDATA_var.ChannelIDs(1)));
for i = 1:length(OFFDATA_var.ChannelIDs)
     OFFtime=[OFFtime; hours(seconds(length(OFFDATA.(OFFDATA_var.ChannelIDs(i)).(string(currentFields(end))))...
         /OFFDATA_var.PNEfs))];
end




% figure
% scatter(OFFDATA.(channels(1)).AllOFFtimes,...
%      repmat(1,1,length(OFFDATA.(channels(1)).AllOFFtimes)))
%  hold on
%  scatter(OFFDATA.(channels(2)).AllOFFtimes,...
%      repmat(2,1,length(OFFDATA.(channels(2)).AllOFFtimes)))
%  scatter(OFFDATA.(channels(3)).AllOFFtimes,...
%      repmat(3,1,length(OFFDATA.(channels(3)).AllOFFtimes)))
%  scatter(OFFDATA.(channels(4)).AllOFFtimes,...
%      repmat(4,1,length(OFFDATA.(channels(4)).AllOFFtimes)))
%  
% binEdge=[0:2:200 205:5:250 260:10:350]; 
% figure
% subplot(2,2,1)
% histogram(((diff(OFFDATA.(channels(1)).nr,1,2)+1/498.2462))*1000,binEdge)
% 
% subplot(2,2,2)
% histogram(((diff(OFFDATA.(channels(2)).nr,1,2)+1/498.2462))*1000,binEdge)
% 
% subplot(2,2,3)
% histogram(((diff(OFFDATA.(channels(3)).nr,1,2)+1/498.2462))*1000,binEdge)
% 
% subplot(2,2,4)
% histogram(((diff(OFFDATA.(channels(4)).nr,1,2)+1/498.2462))*1000,binEdge)
% 
% figure
% v=violinplot([(diff(OFFDATA.(channels(1)).nr,1,2)+1/498.2462)*1000;...
%          (diff(OFFDATA.(channels(2)).nr,1,2)+1/498.2462)*1000;...
%          (diff(OFFDATA.(channels(3)).nr,1,2)+1/498.2462)*1000;...
%          (diff(OFFDATA.(channels(4)).nr,1,2)+1/498.2462)*1000],...
%         [repmat(channelNums{1},length(OFFDATA.(channels(1)).nr),1);...
%          repmat(channelNums{2},length(OFFDATA.(channels(2)).nr),1);...
%          repmat(channelNums{3},length(OFFDATA.(channels(3)).nr),1);...
%          repmat(channelNums{4},length(OFFDATA.(channels(4)).nr),1)],...
%          'ShowData',false);
%      
% figure
% v2=violinplot([coherenceMat(:,1);...
%          coherenceMat(:,2);...
%          coherenceMat(:,3);...
%          coherenceMat(:,4)],...
%         [repmat(channelNums{1},length(coherenceMat),1);...
%          repmat(channelNums{2},length(coherenceMat),1);...
%          repmat(channelNums{3},length(coherenceMat),1);...
%          repmat(channelNums{4},length(coherenceMat),1)],...
%          'ShowData',false);
%      
%  
% 
% [median(((diff(OFFDATA.(channels(1)).nr,1,2)+1/498.2462))*1000),...
%     median(((diff(OFFDATA.(channels(2)).nr,1,2)+1/498.2462))*1000)...
%     median(((diff(OFFDATA.(channels(3)).nr,1,2)+1/498.2462))*1000)...
%     median(((diff(OFFDATA.(channels(4)).nr,1,2)+1/498.2462))*1000)]
% 
% [mean(((diff(OFFDATA.(channels(1)).nr,1,2)+1/498.2462))*1000),...
%     mean(((diff(OFFDATA.(channels(2)).nr,1,2)+1/498.2462))*1000)...
%     mean(((diff(OFFDATA.(channels(3)).nr,1,2)+1/498.2462))*1000)...
%     mean(((diff(OFFDATA.(channels(4)).nr,1,2)+1/498.2462))*1000)]

plotBG = uibuttongroup(g.Channelstats,...
    'Position',[0.75 0.2 0.22 0.6 ],'Visible','off',...
    'SelectionChangedFcn',@CHANNELSTAT_PLOT);
uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Off period durations',...
                  'Units','normalized',...
                  'Position',[0.1 0.8 1 0.2],...
                  'FontSize',15,...
                  'Tag','1',...
                  'HandleVisibility','off');
uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Off period number',...
                  'Units','normalized',...
                  'Position',[0.1 0.6 1 0.2],...
                  'FontSize',15,...
                  'Tag','2',...
                  'HandleVisibility','off');
uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Off period occupancy',...
                  'Units','normalized',...
                  'Position',[0.1 0.4 1 0.2],...
                  'FontSize',15,...
                  'Tag','3',...
                  'HandleVisibility','off');

uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Channel coherence',...
                  'Units','normalized',...
                  'Position',[0.1 0.2 1 0.2],...
                  'FontSize',15,...
                  'Tag','4',...
                  'HandleVisibility','off');              
plotBG.Visible='on';

function CHANNELSTAT_PLOT(source,event)
       cla(findobj('parent',gcbf,'type', 'Axes'),'reset') 
        
       if event.NewValue.Tag=='1'
            violinplot(OFFdurations,OFFdurationsID,...
            'ShowData',false,'ViolinAlpha',0.3);
            ylabel('Off period duration (ms)');
            set(findobj('parent',gcf,'type', 'Axes'),'YScale','log')
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))

       elseif event.NewValue.Tag=='2'
            scatter(OFFDATA_var.ChannelNumsOnly,OFFpnum,[],[0.5 0.5 0.5]);
            set(gca, 'XTick', 1:length(OFFDATA_var.ChannelNumsOnly));
            xlim([0 length(OFFDATA_var.ChannelNumsOnly)+1])
            ylabel('Off period number');
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))

       elseif event.NewValue.Tag=='3'
            scatter(OFFDATA_var.ChannelNumsOnly,OFFtime,[],[0.5 0.5 0.5]);
            set(gca, 'XTick', 1:length(OFFDATA_var.ChannelNumsOnly));
            xlim([0 length(OFFDATA_var.ChannelNumsOnly)+1])
            ylabel('Off period occupancy (hours)');
            
       elseif event.NewValue.Tag=='4'
            violinplot(coherence,coherenceID,...
            'ShowData',false,'ViolinAlpha',0.3);
             ylabel('Channel coherence score');
       end

        
end
end