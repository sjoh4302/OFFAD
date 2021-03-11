function [OFFDATA]=OFFAD_channelstats(OFFDATA,plotType)
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

%% Decide if plotting original or adjusted OFF period data
if plotType==1;
    if isfield(OFFDATA,'StartOPadjusted')
        selectData='OPadjusted';  
    else
        selectData='OP';
    end
else
    selectData='OP';
end

%%% Display plotting settings
if contains(selectData,'adjusted')
    displayName='Adjusted';
else
    displayName='Original';
end
uicontrol(g.Channelstats,'Style', 'text','String',displayName,...
    'FontSize',20,...
    'Units','normalized',...
    'Position',[0.72 0.8 0.28 0.1])
    
%%% Add switch data type button
if plotType==0
    switchName='Switch to adjusted';
    switchCallback='close(findobj(''Tag'',''OFFAD_CHANNELSTATS''));[OFFDATA]=OFFAD_channelstats(OFFDATA,1);';
    visibility='on';
elseif contains(selectData,'adjusted') 
    switchName='Switch to original';
    switchCallback='close(findobj(''Tag'',''OFFAD_CHANNELSTATS''));[OFFDATA]=OFFAD_channelstats(OFFDATA,0);';
    visibility='on';
else
    switchName='Switch to original';
    switchCallback='close(findobj(''Tag'',''OFFAD_CHANNELSTATS''));[OFFDATA]=OFFAD_channelstats(OFFDATA,0);';
    visibility='off';
end

uicontrol(g.Channelstats,'Style', 'pushbutton','String',switchName,...
'FontSize',8,...
'Visible',visibility,...
'Units','normalized',...
'Callback',switchCallback,...
'Position',[0.79 0.8 0.14 0.04])

%% Compute statistics

%%% Generate temporary variable
exampleObject = matfile(OFFDATA.PNEpathin);
PNElength=size(exampleObject,OFFDATA.ChannelsFullName(1),2);
PNEtimeTemp=[1/OFFDATA.PNEfs:1/OFFDATA.PNEfs:PNElength/OFFDATA.PNEfs]';
clear PNElength exampleObject

%%%%%%% Off durations plot
OFFdurations=[];
OFFdurationsID=[];
for i = 1:length(OFFDATA.Channels)
   start_end=[PNEtimeTemp(find(OFFDATA.(['Start',selectData])(:,i)==1)),...
       PNEtimeTemp(find(OFFDATA.(['End',selectData])(:,i)==1))];
   OFFdurations=[OFFdurations; (diff(start_end,1,2)+1/OFFDATA.PNEfs)*1000];
   OFFdurationsID=[OFFdurationsID; repmat(OFFDATA.Channels(i),length(start_end),1)];
   %Store summary info
   OFFDATA.(['Stats',selectData]).MeanDuration(i,1)=mean((diff(start_end,1,2)+1/OFFDATA.PNEfs)*100);
   clear start_end
end

%%%%%%% Channel coherence plot 
for i = 1:length(OFFDATA.Channels)
    for j = 1:length(OFFDATA.Channels)
       coherenceMat(i,j)=length(intersect(PNEtimeTemp(OFFDATA.(['All',selectData])(:,i)),PNEtimeTemp(OFFDATA.(['All',selectData])(:,j))))...
           /sum(OFFDATA.(['All',selectData])(:,i));
    
    end
end
coherenceMat(coherenceMat==1)=NaN;
coherence=reshape(coherenceMat,[],1);
coherenceID=reshape(repmat(OFFDATA.Channels,length(OFFDATA.Channels),1),[],1);
%Store summary info
OFFDATA.(['Stats',selectData]).MeanCoherence=nanmean(coherenceMat,2);
clear coherenceMat

%%%%%%%%% Off period number
OFFDATA.(['Stats',selectData]).OPnumber=full(sum(OFFDATA.(['Start',selectData])))';

%%%%%%%%% Off period occupancy time
OFFDATA.(['Stats',selectData]).OPoccupancy_time=hours(seconds(sum(OFFDATA.(['All',selectData])/OFFDATA.PNEfs)))';

%%%%%%%%% Optional LFP info
try
    %Filter settings
    p1=0.5; p2=100; s1=0.1; s2=120;
    Wp=[p1 p2]/(OFFDATA.LFPfs/2); Ws=[s1 s2]/(OFFDATA.LFPfs/2); Rp=3; Rs=30;
    [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
    [bb1,aa1]=cheby2(n,Rs,Wn);
    LFPamp=[];
    LFPampID=[];
    for i = 1:length(OFFDATA.Channels)
         sig=load(OFFDATA.LFPpathin,OFFDATA.ChannelsFullName(i));
         sig=sig.(OFFDATA.ChannelsFullName(i));
         if OFFDATA.FiltLFP==1
             sig=filtfilt(bb1,aa1,double(sig));
             %DC offset
             sig = sig-mean(sig);
       
         end
         
         try
             LFPampTmp=single(sig(unique(round(mod(PNEtimeTemp(OFFDATA.(['All',selectData])(:,i)),1/OFFDATA.LFPfs)...
                 +PNEtimeTemp(OFFDATA.(['All',selectData])(:,i))/(1/OFFDATA.LFPfs)))))';
             LFPampTmp=LFPampTmp(round(length(LFPampTmp)/2000):end-round(length(LFPampTmp)/2000)); %Ignore edge LFP filtering effects by using middle 99.9% of OFF periods
             LFPamp=[LFPamp;LFPampTmp];
             LFPampID=[LFPampID;repmat(OFFDATA.Channels(i),...
                 length(LFPampTmp),1)];
         end
         
         %Store summary info
         OFFDATA.(['Stats',selectData]).MeanLFPamp(i,1)=mean(sig(unique(round(mod(PNEtimeTemp(OFFDATA.(['All',selectData])(:,i)),1/OFFDATA.LFPfs)...
             +PNEtimeTemp(OFFDATA.(['All',selectData])(:,i))/(1/OFFDATA.LFPfs)))));
        clear sig LFPampTmp
    end
end    
    

%Assess outliers
allChannelstats=[];
channelStatsFields=fields(OFFDATA.(['Stats',selectData]));
channelStatsFields=string(channelStatsFields);
for i = 1:length(channelStatsFields)
    allChannelstats(:,i)=OFFDATA.(['Stats',selectData]).(channelStatsFields(i));
end
%Find nan channels
useChan=find(~isnan(allChannelstats(:,1)));

if size(allChannelstats,2)<length(useChan)
    for i = 1:size(allChannelstats,2)
        mahalStore(:,i)=mahal(allChannelstats(useChan,i),allChannelstats(useChan,i));
    end
    OFFDATA.(['Stats',selectData]).MahalDist=nan(1,size(allChannelstats,1));
    OFFDATA.(['Stats',selectData]).MahalDist(useChan)=mean(mahalStore,2);
    %OFFDATA.Stats.MahalDist=mahal(allChannelstats(:,1:5),allChannelstats(:,1:5));
else
    OFFDATA.(['Stats',selectData]).MahalDist=ones(size(allChannelstats,1),1);
end
%end    

%% Draw figure
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

 if isempty(OFFDATA.LFPpathin)==0
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
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.(['Stats',selectData]).MahalDist,...
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
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.(['Stats',selectData]).MahalDist,...
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
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.(['Stats',selectData]).OPnumber,...
                 [1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            set(gca, 'XTick', 1:length(OFFDATA.Channels));
            xlim([0 length(OFFDATA.Channels)+1])
            set(gca, 'XTickLabels', OFFDATA.Channels)
            ylabel('Off period number');
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))

       elseif event.NewValue.Tag=='4'
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.(['Stats',selectData]).OPoccupancy_time,...
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
            %gscatter([1:length(OFFDATA.Channels)],OFFDATA.Stats.MeanLFPamp,...
            %[1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            %set(gca, 'XTick', 1:length(OFFDATA.Channels));
            %xlim([0 length(OFFDATA.Channels)+1])
            %set(gca, 'XTickLabels', OFFDATA.Channels)
            %ylabel('LFP amplitude (uV)')
            violinplot(LFPamp,LFPampID,...
            'ShowData',false,'ViolinAlpha',0.3);
             ylabel('LFP amplitude(uV)');  
             xlim([0 length(OFFDATA.Channels)+1])
       end

        
end

end