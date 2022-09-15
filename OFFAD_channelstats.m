function [OFFDATA]=OFFAD_channelstats(OFFDATA,plotType)
%
% Channel statistics page: Compute OFF/ON period statistics for each channel
%
% Author: Christian Harding 2022
% OFF Period Automated Detection (OFFAD) toolbox
% christian.harding@sjc.ox.uk
%
% Requires:
% - Statistics and Machine learning toolbox
% - Signal processing toolbox
%


% Get screensize to set figure position
screensize = get( groot, 'Screensize' ); 
    
% Create figure for scrolling window 
g.Channelstats = figure('name', 'OFFAD (OFF_period Automated Detection) - Channel Statistics', ... 
	'numbertitle', 'off', ...
	'Position',[screensize(3)/6,screensize(4)/8,4*screensize(3)/6,screensize(4)/1.4], ...
    'Tag','OFFAD_CHANNELSTATS');
drawnow 

% Plotting original or adjusted (with miniumum duration/maxiumum interruption
% thresholds) OFF period data
if plotType==1;
    if isfield(OFFDATA,'StartOPadjusted')
        selectData='OPadjusted';  
    else
        selectData='OP';
    end
else
    selectData='OP';
end

% Display current plotting selection
if contains(selectData,'adjusted')
    displayName='Adjusted';
else
    displayName='Original';
end
uicontrol(g.Channelstats,'Style', 'text','String',displayName,...
    'FontSize',20,...
    'Units','normalized',...
    'Position',[0.72 0.8 0.28 0.1])
    
% Button to switch between original and adjusted data types
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

% Generate temporary variable to establish MUA length
exampleObject = matfile(OFFDATA.MUApathin);
MUAlength=size(exampleObject,OFFDATA.ChannelsFullName(1),2);
MUAtimeTemp=[1/OFFDATA.MUAfs:1/OFFDATA.MUAfs:MUAlength/OFFDATA.MUAfs]';
clear MUAlength exampleObject

% Statistics 1: OFF period duration
% Find duration of all OFF periods in recording for each channel (ms)
OFFdurations=[];
OFFdurationsID=[];
for i = 1:length(OFFDATA.Channels)
   start_end=[MUAtimeTemp(find(OFFDATA.(['Start',selectData])(:,i)==1)),...
       MUAtimeTemp(find(OFFDATA.(['End',selectData])(:,i)==1))];
   OFFdurations=[OFFdurations; (diff(start_end,1,2)+1/OFFDATA.MUAfs)*1000];
   OFFdurationsID=[OFFdurationsID; repmat(OFFDATA.Channels(i),length(start_end),1)];
   %Store summary info
   OFFDATA.(['Stats',selectData]).MeanDuration(i,1)=mean((diff(start_end,1,2)+1/OFFDATA.MUAfs)*100);
   clear start_end
end

% Statistics 2: OFF period inter-channel coherence
% Calculate average proportion of OFF periods temporally shared between 
% each pair of channels
for i = 1:length(OFFDATA.Channels)
    display(['Calculating channel ',num2str(OFFDATA.Channels(i))])
    for j = 1:length(OFFDATA.Channels)
       coherenceMat(i,j)=(length(intersect(MUAtimeTemp(OFFDATA.(['All',selectData])(:,i)),MUAtimeTemp(OFFDATA.(['All',selectData])(:,j))))...
           /sum(OFFDATA.(['All',selectData])(:,i))...
           +length(intersect(MUAtimeTemp(OFFDATA.(['All',selectData])(:,j)),MUAtimeTemp(OFFDATA.(['All',selectData])(:,i))))...
           /sum(OFFDATA.(['All',selectData])(:,j)))/2;
     end    
end
coherenceMat(coherenceMat==1)=NaN;
coherence=reshape(coherenceMat,[],1);
coherenceID=reshape(repmat(OFFDATA.Channels,length(OFFDATA.Channels),1),[],1);
% Store summary info
OFFDATA.(['Stats',selectData]).MeanCoherence=nanmean(coherenceMat,2);
clear coherenceMat

% Statistics 3: OFF period number
% Find total number of OFF periods in recording for each channel 
OFFDATA.(['Stats',selectData]).OPnumber=full(sum(OFFDATA.(['Start',selectData])))';

% Statistics 4: OFF period occupancy
% Find total time of spent in OFF state in recording for each channel (hours)
OFFDATA.(['Stats',selectData]).OPoccupancy_time=hours(seconds(sum(OFFDATA.(['All',selectData])/OFFDATA.MUAfs)))';

% Statistics 5: LFP profile [OPTIONAL]
% Find average amplitude of LFP associated with each OFF period per channel
try
    % Filter LFP
    p1=0.5; p2=100; s1=0.1; s2=120;
    Wp=[p1 p2]/(OFFDATA.LFPfs/2); Ws=[s1 s2]/(OFFDATA.LFPfs/2); Rp=3; Rs=30;
    [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
    [bb1,aa1]=cheby2(n,Rs,Wn);
    
    % Create empty variables to store LFP ampltiude
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
         
         % Extract concommitant LFP during each MUA OFF period
         try
             LFPampTmp=single(sig(unique(round(mod(MUAtimeTemp(OFFDATA.(['All',selectData])(:,i)),1/OFFDATA.LFPfs)...
                 +MUAtimeTemp(OFFDATA.(['All',selectData])(:,i))/(1/OFFDATA.LFPfs)))))';
             LFPampTmp=LFPampTmp(round(length(LFPampTmp)/2000):end-round(length(LFPampTmp)/2000)); %Ignore edge LFP filtering effects by using middle 99.9% of OFF periods
             LFPamp=[LFPamp;LFPampTmp];
             LFPampID=[LFPampID;repmat(OFFDATA.Channels(i),...
                 length(LFPampTmp),1)];
         end
         
         %Store summary info
         OFFDATA.(['Stats',selectData]).MeanLFPamp(i,1)=mean(sig(unique(round(mod(MUAtimeTemp(OFFDATA.(['All',selectData])(:,i)),1/OFFDATA.LFPfs)...
             +MUAtimeTemp(OFFDATA.(['All',selectData])(:,i))/(1/OFFDATA.LFPfs)))));
        clear sig LFPampTmp
    end
end    
    

% Store summary information
allChannelstats=[];
channelStatsFields=fields(OFFDATA.(['Stats',selectData]));
channelStatsFields=string(channelStatsFields);
for i = 1:length(channelStatsFields)
    allChannelstats(:,i)=OFFDATA.(['Stats',selectData]).(channelStatsFields(i));
end
% Find empty channels to remove
useChan=find(~isnan(allChannelstats(:,1)));


% Statistics 6: Mahalanobis distance
% Calculate Mahalanobis distance for each channels as a measure of
% divergence to help identify outlier channels (e.g. poor clustering)
if length(useChan)>2
    for i = 1:size(allChannelstats,2)-1
        mahalStore(:,i)=mahal(allChannelstats(useChan,i),allChannelstats(useChan,i));
    end
    OFFDATA.(['Stats',selectData]).MahalDist=nan(1,size(allChannelstats,1));
    OFFDATA.(['Stats',selectData]).MahalDist(useChan)=mean(mahalStore,2);
else
    OFFDATA.(['Stats',selectData]).MahalDist=ones(size(allChannelstats,1),1);
end
    

% Buttons to select channels statistic to view
plotBG = uibuttongroup(g.Channelstats,...
    'Position',[0.75 0.2 0.22 0.6 ],'Visible','off',...
    'SelectionChangedFcn',@CHANNELSTAT_PLOT);

uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','OFF period duration',...
                  'Units','normalized',...
                  'Position',[0.1 0.75 1 0.2],...
                  'FontSize',13,...
                  'Tag','1',...
                  'HandleVisibility','off');
uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','OFF period number',...
                  'Units','normalized',...
                  'Position',[0.1 0.6 1 0.2],...
                  'FontSize',13,...
                  'Tag','2',...
                  'HandleVisibility','off');
uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','OFF period occupancy',...
                  'Units','normalized',...
                  'Position',[0.1 0.45 1 0.2],...
                  'FontSize',13,...
                  'Tag','3',...
                  'HandleVisibility','off');

uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Channel coherence',...
                  'Units','normalized',...
                  'Position',[0.1 0.3 1 0.2],...
                  'FontSize',13,...
                  'Tag','4',...
                  'HandleVisibility','off');   
              
uicontrol(plotBG,'Style',...
                  'radiobutton',...
                  'String','Mahalanobis distance',...
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

% Clear temporary variables
clear MUAtimeTemp


% Plot default graph (duration)
plottingColors=[0    0.4470    0.7410
                       0.3010    0.7450    0.9330
                       0.4940    0.1840    0.5560
                       0.8500    0.3250    0.0980
                       0.6350    0.0780    0.1840
                       0.4660    0.6740    0.1880
                       0.9290    0.6940    0.1250];
subplot('position',[0.07 0.12 0.65 0.8])
set(gca,'ColorOrder',plottingColors)
violinplot(OFFdurations,OFFdurationsID,...
            'ShowData',false,'ViolinAlpha',0.3);
            ylabel('Off period duration (ms)');
            set(findobj('parent',gcf,'type', 'Axes'),'YScale','log')
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))
            xlim([0 length(OFFDATA.Channels)+1])

    
% Exit button to return to main output page
uicontrol(g.Channelstats,'Style', 'pushbutton','String','Done',...
    'FontWeight','bold','FontSize',12,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Position',[0.75 0.05 0.2 0.1],'Callback',...
    'set(findobj(''tag'',''OFFAD''),''Visible'',''on'');close(findobj(''tag'',''OFFAD_CHANNELSTATS''));');

    
% Function to change statistic being viewed
function CHANNELSTAT_PLOT(source,event)
       cla(findobj('parent',gcbf,'type', 'Axes'),'reset') 
       
       % Plot duration
       if event.NewValue.Tag=='1'
            violinplot(OFFdurations,OFFdurationsID,...
            'ShowData',false,'ViolinAlpha',0.3);
            ylabel('Off period duration (ms)');
            set(findobj('parent',gcf,'type', 'Axes'),'YScale','log')
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))
            xlim([0 length(OFFDATA.Channels)+1])
       
       % Plot number
       elseif event.NewValue.Tag=='2'
            set(gca,'ColorOrder',plottingColors)
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.(['Stats',selectData]).OPnumber,...
                 [1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            set(gca, 'XTick', 1:length(OFFDATA.Channels));
            xlim([0 length(OFFDATA.Channels)+1])
            set(gca, 'XTickLabels', OFFDATA.Channels)
            ylabel('Off period number');
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))
        
       % Plot occupancy time
       elseif event.NewValue.Tag=='3'
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.(['Stats',selectData]).OPoccupancy_time,...
                [1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            set(gca, 'XTick', 1:length(OFFDATA.Channels));
            xlim([0 length(OFFDATA.Channels)+1])
            set(gca, 'XTickLabels', OFFDATA.Channels)
            ylabel('Off period occupancy (hours)');
        
       % Plot coherence
       elseif event.NewValue.Tag=='4'
            violinplot(coherence,coherenceID,...
            'ShowData',false,'ViolinAlpha',0.3);
             ylabel('Channel coherence score');
             xlim([0 length(OFFDATA.Channels)+1])
       
       % Plot Mahalanobis distance
       elseif   event.NewValue.Tag=='5'
            set(gca,'ColorOrder',plottingColors)
            gscatter([1:length(OFFDATA.Channels)],OFFDATA.(['Stats',selectData]).MahalDist,...
                 [1:length(OFFDATA.Channels)],plottingColors,[],30,'off');
            set(gca, 'XTick', 1:length(OFFDATA.Channels));
            xlim([0 length(OFFDATA.Channels)+1])
            set(gca, 'XTickLabels', OFFDATA.Channels)
            ylabel('Mahalanobis distance');
            set(findobj('parent',gcf,'type', 'Axes'),'YTickLabel',get(findobj('parent',gcf,'type', 'Axes'),'YTick'))
       
       % Plot LFP amplitude
       elseif event.NewValue.Tag=='6'
             violinplot(LFPamp,LFPampID,...
            'ShowData',false,'ViolinAlpha',0.3);
             ylabel('LFP amplitude(uV)');  
             xlim([0 length(OFFDATA.Channels)+1])
       end
       
end

end