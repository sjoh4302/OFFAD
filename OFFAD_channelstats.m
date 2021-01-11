function [OFFDATA]=OFFAD_channelstats(OFFDATA)

%Close loading window
g.Cluster = findobj('tag', 'OFFAD_CLUSTER');
close(g.Cluster)

%Plotting
g.Channelstats = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576/2 432], ...
	'PaperUnits','points', ...
	'name', 'OFFAD (OFF_period Automated Detection) - Channel Statistics', ... 
	'numbertitle', 'off', ...
	'Position',[200 100 800 450], ...
    'Tag','OFFAD_CHANNELSTATS');
drawnow 

%%% Generate temporary variable
exampleObject = matfile(OFFDATA.PNEpathin);
PNElength=size(exampleObject,OFFDATA.ChannelsFullName(1),2);
PNEtimeTemp=[1/OFFDATA.PNEfs:1/OFFDATA.PNEfs:PNElength/OFFDATA.PNEfs]';
clear PNElength exampleObject

%%%%%%% Off durations plot
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
OFFDATA.nr.OPnumber=full(sum(OFFDATA.nr.StartOP))';

%%%%%%%%% Off period occupancy time
OFFDATA.nr.OPoccupancy_time=hours(seconds(sum(OFFDATA.nr.AllOP/OFFDATA.PNEfs)))';


%%%% Make button selection
plotBG = uibuttongroup(g.Channelstats,...
    'Position',[0.75 0.2 0.22 0.6 ],'Visible','off',...
    'SelectionChangedFcn',@CHANNELSTAT_PLOT);

uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Mahalanobis distance',...
                  'Units','normalized',...
                  'Position',[0.1 0.75 1 0.2],...
                  'FontSize',13,...
                  'Tag','1',...
                  'HandleVisibility','off');
uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Off period durations',...
                  'Units','normalized',...
                  'Position',[0.1 0.6 1 0.2],...
                  'FontSize',13,...
                  'Tag','2',...
                  'HandleVisibility','off');
uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Off period number',...
                  'Units','normalized',...
                  'Position',[0.1 0.45 1 0.2],...
                  'FontSize',13,...
                  'Tag','3',...
                  'HandleVisibility','off');
uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Off period occupancy',...
                  'Units','normalized',...
                  'Position',[0.1 0.3 1 0.2],...
                  'FontSize',13,...
                  'Tag','4',...
                  'HandleVisibility','off');

uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Channel coherence',...
                  'Units','normalized',...
                  'Position',[0.1 0.15 1 0.2],...
                  'FontSize',13,...
                  'Tag','5',...
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
                  'FontSize',13,...
                  'Tag','6',...
                  'HandleVisibility','off'); 
              
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
if size(allChannelstats,2)<size(allChannelstats,1)
    OFFDATA.nr.MahalDist=mahal(allChannelstats,allChannelstats);
else
    OFFDATA.nr.MahalDist=ones(size(allChannelstats,1),1);
end
    
%%%% Plot default graph
plottingColors=[0    0.4470    0.7410
                       0.3010    0.7450    0.9330
                       0.4940    0.1840    0.5560
                       0.8500    0.3250    0.0980
                       0.6350    0.0780    0.1840
                       0.4660    0.6740    0.1880
                       0.9290    0.6940    0.1250];
subplot('position',[0.07 0.12 0.65 0.8])
set(gca,'ColorOrder',plottingColors)
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.nr.MahalDist,...
                 [1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            set(gca, 'XTick', 1:length(OFFDATA.Channels));
            xlim([0 length(OFFDATA.Channels)+1])
            set(gca, 'XTickLabels', OFFDATA.Channels)
            ylabel('Mahalanobis distance');
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))
            

    
%%% Exit button
uicontrol(g.Channelstats,'Style', 'pushbutton','String','Done',...
    'FontWeight','bold','FontSize',12,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Position',[0.75 0.05 0.2 0.1],'Callback',...
    'set(findobj(''tag'',''OFFAD''),''Visible'',''on'');close(findobj(''tag'',''OFFAD_CHANNELSTATS''));');

    
function CHANNELSTAT_PLOT(source,event)
       cla(findobj('parent',gcbf,'type', 'Axes'),'reset') 
       
       if   event.NewValue.Tag=='1'
            set(gca,'ColorOrder',plottingColors)
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.nr.MahalDist,...
                 [1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            set(gca, 'XTick', 1:length(OFFDATA.Channels));
            xlim([0 length(OFFDATA.Channels)+1])
            set(gca, 'XTickLabels', OFFDATA.Channels)
            ylabel('Mahalanobis distance');
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))
            
       elseif event.NewValue.Tag=='2'
            violinplot(OFFdurations,OFFdurationsID,...
            'ShowData',false,'ViolinAlpha',0.3);
            ylabel('Off period duration (ms)');
            set(findobj('parent',gcf,'type', 'Axes'),'YScale','log')
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))
            xlim([0 length(OFFDATA.Channels)+1])
            
       elseif event.NewValue.Tag=='3'
            set(gca,'ColorOrder',plottingColors)
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.nr.OPnumber,...
                 [1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            set(gca, 'XTick', 1:length(OFFDATA.Channels));
            xlim([0 length(OFFDATA.Channels)+1])
            set(gca, 'XTickLabels', OFFDATA.Channels)
            ylabel('Off period number');
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))

       elseif event.NewValue.Tag=='4'
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.nr.OPoccupancy_time,...
                [1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            set(gca, 'XTick', 1:length(OFFDATA.Channels));
            xlim([0 length(OFFDATA.Channels)+1])
            set(gca, 'XTickLabels', OFFDATA.Channels)
            ylabel('Off period occupancy (hours)');
            
       elseif event.NewValue.Tag=='5'
            violinplot(coherence,coherenceID,...
            'ShowData',false,'ViolinAlpha',0.3);
             ylabel('Channel coherence score');
             xlim([0 length(OFFDATA.Channels)+1])
             
        elseif event.NewValue.Tag=='6'
             violinplot(LFPamp,LFPampID,...
            'ShowData',false,'ViolinAlpha',0.3);
             ylabel('LFP amplitude');  
             xlim([0 length(OFFDATA.Channels)+1])
       end

        
end
end