function [OFFDATA]=OFFAD_channelstats(OFFDATA)

%Close loading window
g.Cluster = findobj('tag', 'OFFAD_CLUSTER');
set(g.Cluster,'Visible','off')

%Plotting
g.Channelstats = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576/2 432], ...
	'PaperUnits','points', ...
	'name', 'OFFAD (OFF_period Automated Detection)', ... 
	'numbertitle', 'off', ...
	'Position',[200 100 800 450], ...
    'Tag','OFFAD_CHANNELSTATS');

%%% Generate temporary variable
exampleObject = matfile(OFFDATA.PNEpathin);
PNElength=size(exampleObject,OFFDATA.ChannelsFullName(1),2);
PNEtimeTemp=[1/OFFDATA.PNEfs:1/OFFDATA.PNEfs:PNElength/OFFDATA.PNEfs]';
clear PNElength exampleObject

%%%%%%% Off durations plot (default)
OFFdurations=[];
OFFdurationsID=[];
for i = 1:length(OFFDATA.Channels)
   start_end=[PNEtimeTemp(find(OFFDATA.nr.StartOP(:,i)==1)),...
       PNEtimeTemp(find(OFFDATA.nr.EndOP(:,i)==1))];
   OFFdurations=[OFFdurations; (diff(start_end,1,2)+1/OFFDATA.PNEfs)*1000];
   OFFdurationsID=[OFFdurationsID; repmat(OFFDATA.Channels(i),length(start_end),1)];
   %Store summary info
   OFFDATA.nr.MeanDuration(i,1)=mean((diff(start_end,1,2)+1/OFFDATA.PNEfs)*100);
   clear start_end
end

subplot('position',[0.07 0.12 0.65 0.8])
violinplot(OFFdurations,OFFdurationsID,...
'ShowData',false,'ViolinAlpha',0.3);
ylabel('Off period duration (ms)');
set(findobj('parent',gcf,'type', 'Axes'),'YScale','log')
set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))
 

%%%%%%% Channel coherence plot 
for i = 1:length(OFFDATA.Channels)
    for j = 1:length(OFFDATA.Channels)
       coherenceMat(i,j)=length(intersect(PNEtimeTemp(OFFDATA.nr.AllOP(:,i)),PNEtimeTemp(OFFDATA.nr.AllOP(:,j))))...
           /sum(OFFDATA.nr.AllOP(:,i));
    
    end
end
coherenceMat(coherenceMat==1)=NaN;
coherence=reshape(coherenceMat,[],1);
coherenceID=reshape(repmat(OFFDATA.Channels,length(OFFDATA.Channels),1),[],1);
%Store summary info
OFFDATA.nr.MeanCoherence=nanmean(coherenceMat,2);
clear coherenceMat

%%%%%%%%% Off period number
OFFpnum=sum(OFFDATA.nr.StartOP);
OFFDATA.nr.OPnumber=full(sum(OFFDATA.nr.StartOP))';

%%%%%%%%% Off period occupancy time
OFFtime=hours(seconds(sum(OFFDATA.nr.AllOP/OFFDATA.PNEfs)));
OFFDATA.nr.OPoccupancy_time=hours(seconds(sum(OFFDATA.nr.AllOP/OFFDATA.PNEfs)))';



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
%LFP info
try
    LFPamp=[];
    LFPampID=[];
    for i = 1:length(OFFDATA.Channels)
         sig=load(OFFDATA.LFPpathin,OFFDATA.ChannelsFullName(i));
         sig=bandpass(sig.(OFFDATA.ChannelsFullName(i)),[0.5 100],256);
         LFPamp=[LFPamp;single(sig(unique(round(mod(PNEtimeTemp(OFFDATA.nr.AllOP(:,i)),1/OFFDATA.LFPfs)...
             +PNEtimeTemp(OFFDATA.nr.AllOP(:,i))/(1/OFFDATA.LFPfs)))))'];
         LFPampID=[LFPampID;repmat(OFFDATA.Channels(i),...
             length(unique(round(mod(PNEtimeTemp(OFFDATA.nr.AllOP(:,i)),1/OFFDATA.LFPfs)...
              +PNEtimeTemp(OFFDATA.nr.AllOP(:,i))/(1/OFFDATA.LFPfs)))),1)];
         
         %Store summary info
         OFFDATA.nr.MeanLFPamp(i,1)=mean(sig(unique(round(mod(PNEtimeTemp(OFFDATA.nr.AllOP(:,i)),1/OFFDATA.LFPfs)...
             +PNEtimeTemp(OFFDATA.nr.AllOP(:,i))/(1/OFFDATA.LFPfs)))));
        clear sig
    end
    
    
    uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','LFP amplitude',...
                  'Units','normalized',...
                  'Position',[0.1 0 1 0.2],...
                  'FontSize',15,...
                  'Tag','5',...
                  'HandleVisibility','off'); 
              displau('hi')
end              
plotBG.Visible='on';

%Clear temporary variables
clear PNEtimeTemp



%Assess outliers
allChannelstats=[];
channelStatsFields=fields(OFFDATA.nr);
channelStatsFields=string(channelStatsFields(4:end));
for i = 1:length(channelStatsFields)
    allChannelstats(:,i)=OFFDATA.nr.(channelStatsFields(i));
end
MahalDist=mahal(allChannelstats,allChannelstats);


function CHANNELSTAT_PLOT(source,event)
       cla(findobj('parent',gcbf,'type', 'Axes'),'reset') 
       plottingColors=[0    0.4470    0.7410
                       0.3010    0.7450    0.9330
                       0.4940    0.1840    0.5560
                       0.8500    0.3250    0.0980
                       0.6350    0.0780    0.1840
                       0.4660    0.6740    0.1880
                       0.9290    0.6940    0.1250];
                  
       if event.NewValue.Tag=='1'
            violinplot(OFFdurations,OFFdurationsID,...
            'ShowData',false,'ViolinAlpha',0.3);
            ylabel('Off period duration (ms)');
            set(findobj('parent',gcf,'type', 'Axes'),'YScale','log')
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))
            xlim([0 length(OFFDATA.Channels)+1])
            
       elseif event.NewValue.Tag=='2'
            set(gca,'ColorOrder',plottingColors)
            gscatter([1:length(OFFDATA.Channels)],OFFpnum,...
                 [1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            set(gca, 'XTick', 1:length(OFFDATA.Channels));
            xlim([0 length(OFFDATA.Channels)+1])
            set(gca, 'XTickLabels', OFFDATA.Channels)
            ylabel('Off period number');
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))

       elseif event.NewValue.Tag=='3'
            gscatter([1:length(OFFDATA.Channels)],OFFtime,...
                [1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            set(gca, 'XTick', 1:length(OFFDATA.Channels));
            xlim([0 length(OFFDATA.Channels)+1])
            set(gca, 'XTickLabels', OFFDATA.Channels)
            ylabel('Off period occupancy (hours)');
            
       elseif event.NewValue.Tag=='4'
            violinplot(coherence,coherenceID,...
            'ShowData',false,'ViolinAlpha',0.3);
             ylabel('Channel coherence score');
             xlim([0 length(OFFDATA.Channels)+1])
             
        elseif event.NewValue.Tag=='5'
             violinplot(LFPamp,LFPampID,...
            'ShowData',false,'ViolinAlpha',0.3);
             ylabel('LFP amplitude');  
             xlim([0 length(OFFDATA.Channels)+1])
       end

        
end
end