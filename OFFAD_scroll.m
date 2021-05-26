function OFFAD_scroll(OFFDATA); 
%profile on

%%% Figure plotting
g.Scroll = figure('Units','points', ...
    'Color','white',...
	'PaperPosition',[18 180 576/2 432], ...
	'PaperUnits','points', ...
	'name', 'OFFAD (OFF_period Automated Detection) - Scroll', ... 
	'numbertitle', 'off', ...
	'Position',[200 100 800 450], ...
    'UserData',1,...
    'Menubar','none',...
    'Tag','OFFAD_SCROLL');

drawnow

subplot('Position',[0.0 0 0.04 1],'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'XColor','None','YColor','None','Color','None')
xlimPNE=get(gca,'XLim');
ylimPNE=get(gca,'YLim');
ht = text(0.2*xlimPNE(1)+0.2*xlimPNE(2),0.5*ylimPNE(1)+0.5*ylimPNE(2),['PNE (', char(OFFDATA.PNEunit),')']);
set(ht,'Rotation',90)
set(ht,'FontSize',10)
set(ht,'Color',[0, 0.4470, 0.7410])
set(ht,'FontWeight','bold')

subplot('Position',[0.86 0 0.01 1],'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'XColor','None','YColor','None','Color','None')
xlimLFP=get(gca,'XLim');
ylimLFP=get(gca,'YLim');
ht = text(0.8*xlimLFP(1)+0.8*xlimLFP(2),0.5*ylimLFP(1)+0.5*ylimLFP(2),['LFP (', char(OFFDATA.LFPunit),')']);
set(ht,'Rotation',90)
set(ht,'FontSize',10)
set(ht,'Color',[0.8500, 0.3250, 0.0980])
set(ht,'FontWeight','bold')

PNEtmp=load(OFFDATA.PNEpathin,OFFDATA.ChannelsFullName{:});
PNEexampleObject = matfile(OFFDATA.PNEpathin);
PNElength=size(PNEexampleObject,OFFDATA.ChannelsFullName(1),2);

if isempty(OFFDATA.LFPpathin)==0
    %Filter settings
    p1=0.5; p2=100; s1=0.1; s2=120;
    Wp=[p1 p2]/(OFFDATA.LFPfs/2); Ws=[s1 s2]/(OFFDATA.LFPfs/2); Rp=3; Rs=30;
    [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
    [bb1,aa1]=cheby2(n,Rs,Wn);
    LFPtmp=load(OFFDATA.LFPpathin,OFFDATA.ChannelsFullName{:});
    if OFFDATA.FiltLFP==1
        for j = 1:length(OFFDATA.ChannelsFullName)
            LFPtmp.(OFFDATA.ChannelsFullName(j))=filtfilt(bb1,aa1,double(LFPtmp.(OFFDATA.ChannelsFullName(j))));
            %DC offset
            %LFPtmp.(OFFDATA.ChannelsFullName(j))= LFPtmp.(OFFDATA.ChannelsFullName(j))-mean( LFPtmp.(OFFDATA.ChannelsFullName(j)));
        end
    end
    LFPexampleObject = matfile(OFFDATA.LFPpathin);
    LFPlength=size(LFPexampleObject,OFFDATA.ChannelsFullName(1),2);
end


%%%% Create dropdown menu to select horiztonal plot (Hypnogram or average LFP)
uicontrol(g.Scroll,'Style', 'popupmenu','String',{'Hypnogram','Average LFP'},...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.02 0.95 0.08 0.04],...
    'Tag','topPlotMenu',...
    'UserData',0,...
    'Callback',@changeTopPlot);




%%% Generate hypnogram
allEpoch=categorical(nan(ceil(length(OFFDATA.StartOP)/OFFDATA.PNEfs/OFFDATA.epochLen),1));
allCol=repmat([0,0,0],length(allEpoch),1);
vigStateNames={'w','w1','mt','r','r3','nr','nr2','ma'};
vigInf=load(OFFDATA.VSpathin,vigStateNames{:});
allEpoch(vigInf.(vigStateNames{6}))='N'; allEpoch(vigInf.(vigStateNames{7}))='N'; %NREM
allEpoch(vigInf.(vigStateNames{4}))='R'; allEpoch(vigInf.(vigStateNames{5}))='R'; %REM
allEpoch(vigInf.(vigStateNames{1}))='W'; allEpoch(vigInf.(vigStateNames{2}))='W'; allEpoch(vigInf.(vigStateNames{3}))='W'; allEpoch(vigInf.(vigStateNames{8}))='W';%Wake

%Give each state a color
allCol(allEpoch=='W',1)=0; allCol(allEpoch=='W',2)=0.447; allCol(allEpoch=='W',3)=0.7410;
allCol(allEpoch=='R',2)=0.7;
allCol(allEpoch=='N',3)=0;

subplot('Position',[0.12 0.84 0.72 0.14]);
scatter(1:ceil(length(categorical(allEpoch))),categorical(allEpoch),10,allCol,'.')
set(gca,'Tag','topPlot')
set(gca,'NextPlot','Add')
set(gca,'xtick',[])
set(gca,'xcolor',[1,1,1])
set(gca,'xlim',[0,ceil(length(OFFDATA.StartOP)/OFFDATA.PNEfs/OFFDATA.epochLen)])
xline(gca,1,'LineWidth',2,'Color',allCol(1,:))

%Change vigilance state names to full
allEpochFull=allEpoch;
allEpochFull(allEpoch=='W')='WAKE';
allEpochFull(allEpoch=='N')='NREM';
allEpochFull(allEpoch=='R')='REM';

uicontrol(g.Scroll,'Style', 'text','String',string(allEpochFull(1)),...
    'FontWeight','bold','FontSize',23,...
    'Units','normalized',...
    'ForegroundColor',allCol(1,:),... 
    'BackgroundColor',[1 1 1],...
    'Tag','VigState',...
    'Position',[0.15 0.79 0.1 0.05]);  %Plot current state

clear vigStateNames



%%% GUI plotting

%Background for optimization panel
uipanel(g.Scroll,'Position',[0.1 0. 0.7 0.065])

%Select PNE scale
uicontrol(g.Scroll,'Style', 'text','String','PNE max',...
    'FontWeight','bold','FontSize',7,...
    'Units','normalized',...
    'Position',[0.005 0.04 0.06 0.03]);
uicontrol(g.Scroll,'Style', 'edit','String',150,...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.005 0.01 0.06 0.03],...
    'Tag','PNEscale',...
    'Callback',@drawOFFP);

%Select LFP scale
if isempty(OFFDATA.LFPpathin)==0
    uicontrol(g.Scroll,'Style', 'text','String','LFP max',...
    'FontWeight','bold','FontSize',7,...
    'Units','normalized',...
    'Position',[0.82 0.04 0.06 0.03]);
    uicontrol(g.Scroll,'Style', 'edit','String',700,...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.82 0.01 0.06 0.03],...
    'Tag','LFPscale',...
    'Callback',@drawOFFP);
    

end


%Create forward and backward scrolling arrows
uicontrol(g.Scroll,'Style', 'pushbutton','String', '<',...
    'FontWeight','bold','FontSize',16,...
    'Units','normalized',...
    'Position',[0.30 0.8 0.04 0.03],...
    'Callback',@shiftScroll);

uicontrol(g.Scroll,'Style', 'pushbutton','String', '>',...
    'FontWeight','bold','FontSize',16,...
    'Units','normalized',...
    'Position',[0.58 0.8 0.04 0.03],...
     'Callback',@shiftScroll);

%Show X-scale
uicontrol(g.Scroll,'Style', 'text','String','0.5 sec',...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.11 0.09 0.04 0.02]); 


%Create sliders for min duration and max interval
try
    startMinDur=OFFDATA.OFFthresh(1);
    startDurThresh=OFFDATA.OFFthresh;
    
catch
    startMinDur=0;
    startDurThresh=repmat(0,1,length(OFFDATA.Channels));
end
uicontrol(g.Scroll,'Style', 'slider','Value',startMinDur,...
    'Min',0,'Max',250,'SliderStep',[1/250,5/250],...
    'Units','normalized',...
    'Position',[0.11 0.01 0.23 0.03],...
    'Tag','minDur',...
    'Callback',@shiftDur);

uicontrol(g.Scroll,'Style', 'edit','String',startMinDur,...
    'FontSize',8,...
    'Units','normalized',...
    'Position',[0.34 0.01 0.03 0.03],...
    'Tag','minDurVal',...
    'UserData',startDurThresh,...
    'Callback',@shiftDurEdit);

uicontrol(g.Scroll,'Style', 'text','String','Minimum off-period duration (ms)',...
    'FontSize',8,...
    'Units','normalized',...
    'Position',[0.11 0.04 0.26 0.02])

try
    startMaxInt=OFFDATA.ONthresh(1);
    startIntThresh=OFFDATA.ONthresh;
catch
    startMaxInt=0;
    startIntThresh=repmat(0,1,length(OFFDATA.Channels));
end
uicontrol(g.Scroll,'Style', 'slider','Value',startMaxInt,...
    'Min',0,'Max',250,'SliderStep',[1/250,5/250],...
    'Units','normalized',...
    'Position',[0.55 0.01 0.21 0.03],...
    'Tag','maxInter',...
    'Callback',@shiftInter);

uicontrol(g.Scroll,'Style', 'edit','String',startMaxInt,...
    'FontSize',8,...
    'Units','normalized',...
    'Position',[0.76 0.01 0.03 0.03],...
    'Tag','maxInterVal',...
    'UserData',startIntThresh,...
    'Callback',@shiftInterEdit);

uicontrol(g.Scroll,'Style', 'text','String','Maximum off-period interuption (ms)',...
    'FontSize',8,...
    'Units','normalized',...
    'Position',[0.55 0.04 0.21 0.02])


% Create exit buttons
uicontrol(g.Scroll,'Style','pushbutton','String','Store and Exit',...
    'FontWeight','bold','FontSize',12,...
    'Units','normalized',...
    'BackgroundColor',[0.3 0.8 0.8],...
    'Position',[0.89 0.05 0.105 0.04],...
    'Callback',['[OFFDATA]=OFFAD_scrollSave(OFFDATA,'...
    'get(findobj(''Tag'',''minDurVal''),''UserData''),'...
    'get(findobj(''Tag'',''maxInterVal''),''UserData''));'...
    'close(findobj(''Tag'',''OFFAD_SCROLL''))'])

uicontrol(g.Scroll,'Style','pushbutton','String','Exit',...
    'FontWeight','bold','FontSize',12,...
    'Units','normalized',...
    'BackgroundColor',[0.3 0.8 0.8],...
    'Position',[0.89 0.005 0.105 0.04],...
    'Callback','close(findobj(''Tag'',''OFFAD_SCROLL''))')


%Drop down menu to select channel to plot histogram
uicontrol(g.Scroll,'Style', 'popupmenu','String',OFFDATA.ChannelsFullName,...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.91 0.95 0.08 0.04],...
    'Tag','histChan',...
    'Callback',@shiftChan);



%Radiobuttons for channel specific optimizations
try
    if length(unique(OFFDATA.OFFthresh))+length(unique(OFFDATA.ONthresh))==2
        masterVal=1;
        rbVal=repmat(1,1,length(OFFDATA.ChannelsFullName));
        rbEnable='Off';
        chanString='All channels';
    else
        masterVal=0;
        rbVal=[1,repmat(0,1,length(OFFDATA.ChannelsFullName)-1)];
        rbEnable='On';
        chanString='Ch1';
    end
catch
    masterVal=1;
    rbVal=repmat(1,1,length(OFFDATA.ChannelsFullName));
    rbEnable='Off';
    chanString='All channels';
end

%Create master button
uicontrol(g.Scroll,'Style','radiobutton','String','All',...
    'Fontsize',10,...
    'Units','normalized',...    
    'BackgroundColor',[1 1 1],...
    'Value',masterVal,...
    'Tag','RBmaster',...
    'UserData',length(OFFDATA.ChannelsFullName),...
    'Callback',@specificOpt,...
    'Position',[0.03 0.79 0.05 0.03]);  

%Create channel buttons
radioHeight1=linspace(0.09,0.79,length(OFFDATA.ChannelsFullName)+1);
radioHeight2=(radioHeight1(2)-radioHeight1(1))-(radioHeight1(2)-radioHeight1(1))/10;
radioHeight2=radioHeight2-0.027;
radioHeight1=flip(radioHeight1(1:end-1));
for m = 1:length(OFFDATA.ChannelsFullName)
    rbName=['RB',char(string((m)))];
    uicontrol(g.Scroll,'Style','radiobutton',...
    'Units','normalized',...    
    'BackgroundColor',[1 1 1],...
    'Value',rbVal(m),...
    'Enable',rbEnable,...
    'UserData',m,...
    'Tag',rbName,...
    'Callback',@changeOpt,...
    'Position',[0.045 radioHeight1(m)+radioHeight2*0.5+0.007 0.015 0.03])
end

%Display name of channels being opitimzed
uicontrol(g.Scroll,'Style', 'text','String','Optimize OFF periods:',...
    'FontWeight','bold','FontSize',8,...
    'Units','normalized',...
    'Position',[0.37 0.04 0.18 0.02])
uicontrol(g.Scroll,'Style', 'text','String',chanString,...
    'FontWeight','bold','FontSize',10,...
    'Units','normalized',...
    'Tag','optChan',...
    'UserData',1,...
    'Position',[0.37 0.01 0.18 0.03])





%Display current time in seconds
uicontrol(g.Scroll,'Style', 'edit','String','0',...
    'FontWeight','bold','FontSize',10,...
    'Units','normalized',...
    'Position',[0.48 0.8 0.06 0.028],'CreateFcn',@drawOFFP,...
    'Tag','scrollTime',...
    'Callback',@drawOFFP);
uicontrol(g.Scroll,'Style', 'text','String','Recording time(s)',...
    'FontSize',10,...
    'BackgroundColor',[1 1 1],...
    'Units','normalized',...
    'Position',[0.38 0.8 0.1 0.028]);















function drawOFFP(~,~)
oldTime=str2num(get(findobj('Tag','scrollTime'),'String'));
newTime=floor(oldTime/OFFDATA.epochLen)*4; %Round to nearest epoch start time
%Check selected value within bounds
if newTime<0
    newTime=num2str(0);  
end
if newTime>round(length(OFFDATA.AllOP)/OFFDATA.PNEfs)-5
    newTime=num2str(((floor(((length(OFFDATA.AllOP)/OFFDATA.PNEfs)/OFFDATA.epochLen))-1)*4));
end
set(findobj('Tag','scrollTime'),'String',newTime); %set new time

%Reset sliders
set(findobj('Tag','minDur'),'Value',str2num(get(findobj('Tag','minDurVal'),'String')));
set(findobj('Tag','maxInter'),'Value',str2num(get(findobj('Tag','maxInterVal'),'String')));

startSec=str2num(get(findobj('Tag','scrollTime'),'String'));
startPNE=round(startSec/(1/OFFDATA.PNEfs))+1;
endPNE=round(startSec/(1/OFFDATA.PNEfs))+ceil(OFFDATA.PNEfs*4);

tmpPNEtime=[startPNE/OFFDATA.PNEfs:1/OFFDATA.PNEfs:endPNE/OFFDATA.PNEfs];
startLFP=round(mod(startPNE/OFFDATA.PNEfs,1/OFFDATA.LFPfs)+startPNE/OFFDATA.PNEfs/(1/OFFDATA.LFPfs));
endLFP=round(mod(endPNE/OFFDATA.PNEfs,1/OFFDATA.LFPfs)+endPNE/OFFDATA.PNEfs/(1/OFFDATA.LFPfs));
tmpLFPtime=[startLFP/OFFDATA.LFPfs:1/OFFDATA.LFPfs:endLFP/OFFDATA.LFPfs];

%Get plotting scales
PNEscale=str2num(get(findobj('Tag','PNEscale'),'String'));
if isempty(OFFDATA.LFPpathin)==0
    LFPscale=str2num(get(findobj('Tag','LFPscale'),'String'));
end


numChan=length(OFFDATA.ChannelsFullName);
subplot2=linspace(0.1,0.8,numChan+1);
subplot4=(subplot2(2)-subplot2(1))-(subplot2(2)-subplot2(1))/10;
subplot2=flip(subplot2(1:end-1));
for i = 1:numChan
    uicontrol(g.Scroll,'Style', 'text','String',OFFDATA.ChannelsFullName(i),...
    'FontWeight','bold','FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.0135 subplot2(i)+subplot4*0.5-0.02 0.03 0.03]);  %Plot channel name

    subplot('Position',[0.08 subplot2(i) 0.76 subplot4]);
    yyaxis left
    stem(tmpPNEtime,PNEtmp.(OFFDATA.ChannelsFullName(i))(startPNE:endPNE),'-','LineWidth',1,'Marker','none','ShowBaselin','off')
    hold on
    
    
    %%%%Select which OFF periods to plot (duration/interval criteria)
    tmpOFFStarts=find(OFFDATA.StartOP(startPNE:endPNE,i));
    tmpOFFEnds=find(OFFDATA.EndOP(startPNE:endPNE,i));
    
      
    if ~isempty(tmpOFFStarts)
        %%Fix edge window problems
        %Add previous OFF period before window
        try
        tmpOFFStarts=[max(find(OFFDATA.StartOP(:,i),sum((startPNE-...
                    find(OFFDATA.StartOP(:,i)))>0)))-startPNE+1;tmpOFFStarts];
            %If window starts in an ON state     
            if tmpOFFStarts(2)<=tmpOFFEnds(1) 
                 tmpOFFEnds=[max(find(OFFDATA.EndOP(:,i),sum((startPNE-...
                        find(OFFDATA.StartOP(:,i)))>0)))-startPNE+1;tmpOFFEnds];
            end 
        end
        try
        %Add next OFF period after window
        tmpOFFEnds=[tmpOFFEnds;max(find(OFFDATA.EndOP(:,i),length(find(OFFDATA.EndOP(:,i)))-...
                     sum((find(OFFDATA.EndOP(:,i)))-endPNE>1)+1))-startPNE+1];
            %If window ends in ON state
            if tmpOFFEnds(end-1)>=tmpOFFStarts(end)
                tmpOFFStarts=[tmpOFFStarts;max(find(OFFDATA.StartOP(:,i),length(find(OFFDATA.EndOP(:,i)))-...
                         sum((find(OFFDATA.EndOP(:,i)))-endPNE>1)+1))-startPNE+1];
            end
        end
        
        %Remove short gaps
        allIntThresh=get(findobj('Tag','maxInterVal'),'UserData');
        chooseON=tmpOFFStarts(2:end)-tmpOFFEnds(1:end-1)-1<(OFFDATA.PNEfs/1000*allIntThresh(i));
        chooseONstart=~[0;chooseON];
        chooseONend=~[chooseON;0];
        tmpOFFStarts=tmpOFFStarts(chooseONstart);
        tmpOFFEnds=tmpOFFEnds(chooseONend);
        
        %Remove short off-periods
        allDurThresh=get(findobj('Tag','minDurVal'),'UserData');
        chooseOFF=tmpOFFEnds-tmpOFFStarts+1>(OFFDATA.PNEfs/1000*allDurThresh(i));
        tmpOFFStarts=tmpOFFStarts(chooseOFF);
        tmpOFFEnds=tmpOFFEnds(chooseOFF);
        
        %Trim start/end times to window length
        tmpOFFStarts(tmpOFFStarts<1)=1;
        tmpOFFEnds(tmpOFFEnds>ceil(OFFDATA.PNEfs*4))=ceil(OFFDATA.PNEfs*4);
        
        tmpOFFall=[];
        for k = 1:length(tmpOFFStarts)
            tmpOFFall=[tmpOFFall,[tmpOFFStarts(k):tmpOFFEnds(k)]];
        end

        stem(tmpPNEtime(tmpOFFall),PNEtmp.(OFFDATA.ChannelsFullName(i))(startPNE-1+tmpOFFall),'-','LineWidth',1,'Marker','none','ShowBaselin','off','Color','red','LineStyle','-')
    end
    ylim([-PNEscale PNEscale])
    clear chooseOFF tmpOFFStarts tmpOFFEnds tmpOFFall
    
    if isempty(OFFDATA.LFPpathin)==0
        yyaxis right
        plot(tmpLFPtime,LFPtmp.(OFFDATA.ChannelsFullName(i))(startLFP:endLFP),'-')
        ylim([-LFPscale LFPscale])
    end
    
    xlim([tmpPNEtime(1) tmpPNEtime(end)])
    set(gca,'Box','off')
    ax=gca;
    ax.XAxis.Exponent=0;
    
    if i<numChan
        set(gca,'Xcolor','none')
    end
    
end

%Set x-axis to show start/end of segment
set(gca,'XLim',[round(tmpPNEtime(1)),round(tmpPNEtime(end))])
   
%%%%%%% Draw ON-OFF period histograms
try
if get(findobj('Tag','OFFAD_SCROLL'),'UserData')==1;
     % Recalculate OFF start and end time
     for i = get(findobj('Tag','histChan'),'Value')
         adjOFFStarts=find(OFFDATA.StartOP(:,i));
         adjOFFEnds=find(OFFDATA.EndOP(:,i));
        
        %Remove short gaps
        chooseON=adjOFFStarts(2:end)-adjOFFEnds(1:end-1)-1<(OFFDATA.PNEfs/1000*str2num(get(findobj('Tag','maxInterVal'),'String')));
        chooseONstart=~[0;chooseON];
        chooseONend=~[chooseON;0];
        adjOFFStarts=adjOFFStarts(chooseONstart);
        adjOFFEnds=adjOFFEnds(chooseONend);
        
        %Remove short off-periods
        chooseOFF=adjOFFEnds-adjOFFStarts+1>(OFFDATA.PNEfs/1000*str2num(get(findobj('Tag','minDurVal'),'String')));
        adjOFFStarts=adjOFFStarts(chooseOFF);
        adjOFFEnds=adjOFFEnds(chooseOFF);
     end 
     
     %%% Plot OFF periods histrogram
    [oldOFF,oldBins]=histcounts((find(OFFDATA.EndOP(:,i))-find(OFFDATA.StartOP(:,i))+1)*1000/OFFDATA.PNEfs,0:1000/OFFDATA.PNEfs:200);
    [newOFF,newBins]=histcounts((adjOFFEnds-adjOFFStarts+1)*1000/OFFDATA.PNEfs,0:1000/OFFDATA.PNEfs:200);
    subplot('Position',[0.91 0.56 0.08 0.38])
    cla
    a1=area(oldBins(1:end-1)+diff(oldBins)/2,oldOFF,'LineStyle','none');
    a1.FaceAlpha = 0.8;
    hold on
    a2=area(newBins(1:end-1)+diff(newBins)/2,newOFF,'LineStyle','none');
    a2.FaceAlpha = 0.8;
    ylabel('Count')
    ax=gca;
    ax.YScale='log';
    legend('Original','Adjusted')
    text(0.2,0.8,'Off periods','Units','normalized')
    
    %%%Plot ON periods histogram
    [oldON,oldBins]=histcounts((find(OFFDATA.StartOP(:,i),length(find(OFFDATA.StartOP(:,i)))-1,'last')-...
        find(OFFDATA.EndOP(:,i),length(find(OFFDATA.EndOP(:,i)))-1)-1)*1000/OFFDATA.PNEfs,0:1000/OFFDATA.PNEfs:200);
    [newON,newBins]=histcounts((adjOFFStarts(2:end)-adjOFFEnds(1:end-1)-1)*1000/OFFDATA.PNEfs,0:1000/OFFDATA.PNEfs:200);
    subplot('Position',[0.91 0.15 0.08 0.38])
    cla
    a1=area(oldBins(1:end-1)+diff(oldBins)/2,oldON,'LineStyle','none');
    a1.FaceAlpha = 0.8;
    hold on
    a2=area(newBins(1:end-1)+diff(newBins)/2,newON,'LineStyle','none');
    a2.FaceAlpha = 0.8;
    ylabel('Count')
    xlabel('Duration (ms)')
    ax=gca;
    ax.YScale='log';
    text(0.2,0.8,'On periods','Units','normalized')
    
    set(findobj('Tag','OFFAD_SCROLL'),'UserData',0)
    clear chooseON chooseONstart chooseONend adjOFFStarts adjOFFEnds
end
 
catch
     subplot('Position',[0.91 0.56 0.08 0.38])
     cla
     subplot('Position',[0.91 0.15 0.08 0.38])
     cla
end 
 

%%%% Change top plot
if get(findobj('Tag','topPlotMenu'),'UserData')==1
axes(findobj('Tag','topPlot'))
cla(gca)
colorbar(gca,'off')
ylabel(gca,'')

%%% 1: Plot hypnogram
if get(findobj('Tag','topPlotMenu'),'Value')==1
scatter(1:ceil(length(categorical(allEpoch))),categorical(allEpoch),10,allCol,'.')
set(gca,'xtick',[])
set(gca,'xcolor',[1,1,1])
set(gca,'xlim',[0,ceil(length(OFFDATA.StartOP)/OFFDATA.PNEfs/OFFDATA.epochLen)])
xline(gca,1,'LineWidth',2,'Color',allCol(1,:))


%%% 2: Plot average slow wave
elseif  get(findobj('Tag','topPlotMenu'),'Value')==2
swChan=get(findobj('Tag','histChan'),'Value'); %select channel to generate slow wave profile
%%% Recalculate OFF start and end time
adjOFFStarts=find(OFFDATA.StartOP(:,swChan));
adjOFFEnds=find(OFFDATA.EndOP(:,swChan));

%Remove short gaps
chooseON=adjOFFStarts(2:end)-adjOFFEnds(1:end-1)-1<(OFFDATA.PNEfs/1000*str2num(get(findobj('Tag','maxInterVal'),'String')));
chooseONstart=~[0;chooseON];
chooseONend=~[chooseON;0];
adjOFFStarts=adjOFFStarts(chooseONstart);
adjOFFEnds=adjOFFEnds(chooseONend);

%Remove short off-periods
chooseOFF=adjOFFEnds-adjOFFStarts+1>(OFFDATA.PNEfs/1000*str2num(get(findobj('Tag','minDurVal'),'String')));
adjOFFStarts=adjOFFStarts(chooseOFF);
adjOFFEnds=adjOFFEnds(chooseOFF);


%Remove off-periods at start/end of recordings
minOFF=ceil(OFFDATA.PNEfs); maxOFF=length(OFFDATA.StartOP)-ceil(OFFDATA.PNEfs);
fullOFFindex=adjOFFStarts>minOFF&adjOFFStarts<maxOFF;
adjOFFStarts=adjOFFStarts(fullOFFindex);
adjOFFEnds=adjOFFEnds(fullOFFindex);

%Get adjusted off-period durations
MSdurOP=(adjOFFEnds-adjOFFStarts+1)*(1/OFFDATA.PNEfs)*1000; %Off period length in MS
LFPdurOP=ceil(((adjOFFEnds-adjOFFStarts+1)*OFFDATA.LFPfs)/OFFDATA.PNEfs)-1;  %Off period length in LFP data points
clear adjOFFEnds

numOFFPsamp=length(adjOFFStarts);
dispRangeSec=1;%range of average LFP in seconds
dispRangeMSec=dispRangeSec*1000;
dispRangeLFP=(floor(OFFDATA.LFPfs*dispRangeSec)+1+mod(floor(OFFDATA.LFPfs*dispRangeSec),2));%Number of LFP data points to display per OFF period
dispRangeLFPhalf=floor(dispRangeLFP/2);

%Average OFF period LFPs in 10th percentiles
[sortOFFlengths,sortOFFpos]=sort(MSdurOP);
for k=1:10
    currIndex=ceil(length(MSdurOP)*((k-1)/10))+1:ceil(length(MSdurOP)*(k/10));
    avPercentileOFFlength(k)=mean(sortOFFlengths(currIndex));
    centeredOFFperiods=nan(length(currIndex),dispRangeLFP);
    %Extract OFF period LFP
    for swSamp = 1:length(currIndex)
        selectOP=sortOFFpos(currIndex(swSamp));
        LFPstart=round((adjOFFStarts(selectOP)/OFFDATA.PNEfs)*OFFDATA.LFPfs);
        centeredOFFperiods(swSamp,:)=[LFPtmp.(OFFDATA.ChannelsFullName(swChan))(LFPstart-dispRangeLFPhalf:LFPstart+dispRangeLFPhalf)];    
    end
    avePercentileOFF(:,k)=mean(centeredOFFperiods);
    clear centeredOFFperiods
end

%Plot average LFPs
cMap=flip(copper);
lineCols=cMap(round(rescale(avPercentileOFFlength)*255)+1,:);%rescale lengths to colormap
lineCols=flip(lineCols);
avePercentileOFF=flip(avePercentileOFF,2);
LFPtime=linspace(-dispRangeMSec/2,dispRangeMSec/2,dispRangeLFP)';
for i=1:length(avPercentileOFFlength)
    plot(LFPtime,avePercentileOFF(:,i),'Color',lineCols(i,:),'LineWidth',3)
    hold on
end
set(gca,'XAxisLocation', 'top')
set(gca,'xtick',[LFPtime(1):100:LFPtime(end)])
set(gca,'xcolor',[0,0,0])
xlim([LFPtime(1),LFPtime(end)])
ylim([min(min((avePercentileOFF))),max(max((avePercentileOFF)))])
ylabel(char(OFFDATA.LFPunit))
colormap(cMap)
c=colorbar('Ticks',0:0.2:1,...
         'TickLabels',round(linspace(avPercentileOFFlength(1),avPercentileOFFlength(end),6)));
c.Label.String = 'Mean duration (ms)';
clear chooseON chooseONstart chooseONend adjOFFStarts avePercentileOFF avPercentileOFFlength centeredOFFperiods fullOFFindex
end
set(findobj('Tag','topPlotMenu'),'UserData',0);
end

%Plot current epoch vig state info
set(findobj('Tag','VigState'),'String',string(allEpochFull((str2num(get(findobj('Tag','scrollTime'),'String'))/OFFDATA.epochLen)+1)))
set(findobj('Tag','VigState'),'ForegroundColor',allCol((str2num(get(findobj('Tag','scrollTime'),'String'))/OFFDATA.epochLen)+1,:))
if get(findobj('Tag','topPlotMenu'),'Value')==1 
set(findobj(get(findobj('Tag','topPlot'),'Children'),'Type','ConstantLine'),...
    'Value',(str2num(get(findobj('Tag','scrollTime'),'String'))/OFFDATA.epochLen)+1)
set(findobj(get(findobj('Tag','topPlot'),'Children'),'Type','ConstantLine'),...
    'Color',allCol((str2num(get(findobj('Tag','scrollTime'),'String'))/OFFDATA.epochLen)+1,:))
end


end





function shiftScroll(source,~)

if source.String=='<'
    newTime=str2num(get(findobj('Tag','scrollTime'),'String'))-4;
    set(findobj('Tag','scrollTime'),'String',num2str(newTime));
elseif source.String=='>'
    newTime=str2num(get(findobj('Tag','scrollTime'),'String'))+4;
    set(findobj('Tag','scrollTime'),'String',num2str(newTime));
end

%Set toplot to replot if hypnogram
if  get(findobj('Tag','topPlotMenu'),'Value')==1
    set(findobj('Tag','topPlotMenu'),'UserData',1)
end

drawOFFP

end

function shiftDur(source,~) 
newVal=get(findobj('Tag','minDur'),'Value');
set(findobj('Tag','minDurVal'),'String',num2str(newVal))

%Set channel tresholds
if get(findobj('Tag','RBmaster'),'Value')==1
    set(findobj('Tag','minDurVal'),'UserData',repmat(newVal,1,get(findobj('Tag','RBmaster'),'UserData'))); 
else
    chanDurThresh=get(findobj('Tag','minDurVal'),'UserData');
    currChan=get(findobj('Tag','optChan'),'UserData');
    chanDurThresh(currChan)=newVal;
    set(findobj('Tag','minDurVal'),'UserData',chanDurThresh);
end

%Set histograms to replot
set(findobj('Tag','OFFAD_SCROLL'),'UserData',1)
%Set top plot to replot if average LFP
if  get(findobj('Tag','topPlotMenu'),'Value')==2
    set(findobj('Tag','topPlotMenu'),'UserData',1)
end

drawOFFP
end

function shiftDurEdit(~,~)
set(findobj('Tag','minDur'),'Value',str2num(get(findobj('Tag','minDurVal'),'String')));
shiftDur
end




function shiftInter(~,~)
newVal=get(findobj('Tag','maxInter'),'Value');
set(findobj('Tag','maxInterVal'),'String',num2str(newVal))
%Set histograms to replot
set(findobj('Tag','OFFAD_SCROLL'),'UserData',1)
%Set toplot to replot if average LFP
if  get(findobj('Tag','topPlotMenu'),'Value')==2
    set(findobj('Tag','topPlotMenu'),'UserData',1)
end

%Set channel tresholds
if get(findobj('Tag','RBmaster'),'Value')==1
    set(findobj('Tag','maxInterVal'),'UserData',repmat(newVal,1,get(findobj('Tag','RBmaster'),'UserData'))); 
else
    chanInterThresh=get(findobj('Tag','maxInterVal'),'UserData');
    currChan=get(findobj('Tag','optChan'),'UserData');
    chanInterThresh(currChan)=newVal;
    set(findobj('Tag','maxInterVal'),'UserData',chanInterThresh);
end

drawOFFP
end

function shiftInterEdit(~,~)
set(findobj('Tag','maxInter'),'Value',str2num(get(findobj('Tag','maxInterVal'),'String')));
shiftInter
end





function shiftChan(source,~)
%Set histograms to replot
set(findobj('Tag','OFFAD_SCROLL'),'UserData',1)
%Set toplot to replot if average LFP
if  get(findobj('Tag','topPlotMenu'),'Value')==2
    set(findobj('Tag','topPlotMenu'),'UserData',1)
end

drawOFFP
end




function specificOpt(source,~)
if source.Value == 0
    
    %Disable radio buttons
    for j = 1:source.UserData
        rbName=['RB',char(string((j)))];
        set(findobj('Tag',rbName),'Enable','on')
        if j>1
            set(findobj('Tag',rbName),'Value',0)
        end
    end
    set(findobj('Tag','optChan'),'String',OFFDATA.ChannelsFullName(get(findobj('Tag','RB1'),'UserData')))
    set(findobj('Tag','optChan'),'UserData',get(findobj('Tag','RB1'),'UserData'))
    
    %Display correct channel histograms
    set(findobj('Tag','histChan'),'Value',1)
    shiftChan
    
elseif source.Value == 1

    %Warning
    answer = questdlg('Optimization thresholds for all channels will be reset. Do you wish to continue?', ...
	'Warning', ...
	'Yes','No','No');
    waitfor(answer)
    % Handle response
    switch answer
        case 'Yes'
            

        %Disable radio button
        for j = 1:source.UserData
            rbName=['RB',char(string((j)))];
            set(findobj('Tag',rbName),'Enable','off')
            set(findobj('Tag',rbName),'Value',1)
        end
        set(findobj('Tag','optChan'),'String','All channels')
        set(findobj('Tag','optChan'),'UserData',0)

        %Reset optimization thresholds
        set(findobj('Tag','minDurVal'),'String',0)
        set(findobj('Tag','maxInterVal'),'String',0)
        set(findobj('Tag','minDurVal'),'UserData',repmat(0,1,get(findobj('Tag','RBmaster'),'UserData'))); 
        set(findobj('Tag','maxInterVal'),'UserData',repmat(0,1,get(findobj('Tag','RBmaster'),'UserData'))); 

        %Display correct channel histograms
        set(findobj('Tag','histChan'),'Value',1)
        shiftChan
        
    case 'No'
        source.Value=0;  
    end
    
end
end



function changeOpt(source,~)
%Set all channels to 0
chanNum=get(findobj('Tag','RBmaster'),'UserData');
for j = 1:chanNum
        rbName=['RB',char(string((j)))];
        set(findobj('Tag',rbName),'Value',0)
end
%Set current channel to 1
source.Value=1;

%Display current channel being optimized
set(findobj('Tag','optChan'),'String',OFFDATA.ChannelsFullName(source.UserData))
set(findobj('Tag','optChan'),'UserData',source.UserData)

%Display correct channel histograms
set(findobj('Tag','histChan'),'Value',source.UserData)

%Display correct optimization thresholds
allIntThresh=get(findobj('Tag','maxInterVal'),'UserData');
set(findobj('Tag','maxInterVal'),'String',num2str(allIntThresh(source.UserData)))
allDurThresh=get(findobj('Tag','minDurVal'),'UserData');
set(findobj('Tag','minDurVal'),'String',num2str(allDurThresh(source.UserData)))

%Display new histograms
shiftChan

end

function changeTopPlot(source,~)
%Set top plot to replot
set(findobj('Tag','topPlotMenu'),'UserData',1);
drawOFFP
end


end

%clear PNEtmp
%clear LFPtmp
%%%% FAST LOADING EXAMPLE SCRIPT
%save('aa.mat','aa','-v7.3')
%exampleObject = matfile('D:\DPhil\Off-period detection\Matlab scripts\aa.mat');
%PNEtmp=exampleObject.aa(1,1:1000);