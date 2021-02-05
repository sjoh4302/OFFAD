function OFFAD_scroll(OFFDATA); 

%Plotting
g.Scroll = figure('Units','points', ...
    'Color','white',...
	'PaperPosition',[18 180 576/2 432], ...
	'PaperUnits','points', ...
	'name', 'OFFAD (OFF_period Automated Detection) - Scroll', ... 
	'numbertitle', 'off', ...
	'Position',[200 100 800 450], ...
    'Tag','OFFAD_SCROLL');

subplot('Position',[0.0 0 0.06 1],'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'XColor','None','YColor','None','Color','None')
xlimPNE=get(gca,'XLim');
ylimPNE=get(gca,'YLim');
ht = text(0.2*xlimPNE(1)+0.2*xlimPNE(2),0.5*ylimPNE(1)+0.5*ylimPNE(2),['PNE (', char(OFFDATA.PNEunit),')']);
set(ht,'Rotation',90)
set(ht,'FontSize',16)

subplot('Position',[0.94 0 0.06 1],'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'XColor','None','YColor','None','Color','None')
xlimLFP=get(gca,'XLim');
ylimLFP=get(gca,'YLim');
ht = text(0.8*xlimLFP(1)+0.8*xlimLFP(2),0.5*ylimLFP(1)+0.5*ylimLFP(2),['LFP (', char(OFFDATA.LFPunit),')']);
set(ht,'Rotation',90)
set(ht,'FontSize',16)

drawnow

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
            LFPtmp.(OFFDATA.ChannelsFullName(j))=filtfilt(bb1,aa1,LFPtmp.(OFFDATA.ChannelsFullName(j)));
        end
    end
    LFPexampleObject = matfile(OFFDATA.LFPpathin);
    LFPlength=size(LFPexampleObject,OFFDATA.ChannelsFullName(1),2);
end

%Select PNE scale
uicontrol(g.Scroll,'Style', 'text','String','PNE max (uV)',...
    'FontWeight','bold','FontSize',10,...
    'Units','normalized',...
    'Position',[0.01 0.04 0.09 0.03]);
uicontrol(g.Scroll,'Style', 'edit','String',150,...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.01 0.01 0.09 0.03],...
    'Tag','PNEscale',...
    'Callback',@drawOFFP);

%Select LFP scale
if isempty(OFFDATA.LFPpathin)==0
    uicontrol(g.Scroll,'Style', 'text','String','LFP max (uV)',...
    'FontWeight','bold','FontSize',10,...
    'Units','normalized',...
    'Position',[0.9 0.04 0.09 0.03]);
    uicontrol(g.Scroll,'Style', 'edit','String',700,...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.9 0.01 0.09 0.03],...
    'Tag','LFPscale',...
    'Callback',@drawOFFP);
    

end


%Create forward and backward scrolling arrows
uicontrol(g.Scroll,'Style', 'pushbutton','String', '<',...
    'FontWeight','bold','FontSize',20,...
    'Units','normalized',...
    'Position',[0.40 0.01 0.04 0.05],...
    'Callback',@shiftScroll);

uicontrol(g.Scroll,'Style', 'pushbutton','String', '>',...
    'FontWeight','bold','FontSize',20,...
    'Units','normalized',...
    'Position',[0.56 0.01 0.04 0.05],...
     'Callback',@shiftScroll);

%Show X-scale
uicontrol(g.Scroll,'Style', 'text','String','0.5 sec',...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.11 0.09 0.04 0.02]); 


%Create sliders for min duration and max interval
uicontrol(g.Scroll,'Style', 'slider','Value',0,...
    'Min',0,'Max',250,'SliderStep',[1/250,5/250],...
    'Units','normalized',...
    'Position',[0.12 0.01 0.23 0.03],...
    'Tag','minDur',...
    'Callback',@shiftDur);
uicontrol(g.Scroll,'Style', 'edit','String',0,...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.35 0.01 0.02 0.03],...
    'Tag','minDurVal',...
    'Callback',@drawOFFP);
uicontrol(g.Scroll,'Style', 'text','String','Minimum off-period duration (ms)',...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.12 0.04 0.26 0.03])
    
uicontrol(g.Scroll,'Style', 'slider','Value',0,...
    'Min',0,'Max',250,'SliderStep',[1/250,5/250],...
    'Units','normalized',...
    'Position',[0.62 0.01 0.23 0.03],...
    'Tag','maxInter',...
    'Callback',@shiftInter);
uicontrol(g.Scroll,'Style', 'edit','String',0,...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.85 0.01 0.02 0.03],...
    'Tag','maxInterVal',...
    'Callback',@drawOFFP);
uicontrol(g.Scroll,'Style', 'text','String','Maximum off-period interuption (ms)',...
    'FontSize',8,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.62 0.04 0.26 0.03])



%Display current time in seconds
uicontrol(g.Scroll,'Style', 'edit','String','0',...
    'FontWeight','bold','FontSize',10,...
    'Units','normalized',...
    'Position',[0.46 0.01 0.08 0.03],'CreateFcn',@drawOFFP,...
    'Tag','scrollTime',...
    'Callback',@drawOFFP);
uicontrol(g.Scroll,'Style', 'text','String','Time(s)',...
    'FontWeight','bold','FontSize',10,...
    'Units','normalized',...
    'Position',[0.46 0.04 0.08 0.03]);















function drawOFFP(~,~)
oldTime=str2num(get(findobj('Tag','scrollTime'),'String'));
%Check selected value within bounds
if oldTime<0
    newTime=num2str(0);
    set(findobj('Tag','scrollTime'),'String',newTime);
end
if oldTime>round(length(OFFDATA.AllOP)/OFFDATA.PNEfs)-5
    newTime=num2str(round(length(OFFDATA.AllOP)/OFFDATA.PNEfs)-5);
    set(findobj('Tag','scrollTime'),'String',newTime);
end

%Reset sliders
set(findobj('Tag','minDur'),'Value',str2num(get(findobj('Tag','minDurVal'),'String')));
set(findobj('Tag','maxInter'),'Value',str2num(get(findobj('Tag','maxInterVal'),'String')));

startSec=str2num(get(findobj('Tag','scrollTime'),'String'));
startPNE=round(startSec/(1/OFFDATA.PNEfs))+1;
endPNE=round(startSec/(1/OFFDATA.PNEfs))+2000;

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
subplot2=linspace(0.1,0.95,numChan+1);
subplot4=(subplot2(2)-subplot2(1))-(subplot2(2)-subplot2(1))/10;
for i = 1:numChan
    uicontrol(g.Scroll,'Style', 'text','String',OFFDATA.ChannelsFullName(i),...
    'FontWeight','bold','FontSize',10,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position',[0.02 subplot2(i)+subplot4*0.45 0.04 0.03]);  %Plot channel name

    subplot('Position',[0.08 subplot2(i) 0.84 subplot4]);
    yyaxis left
    stem(tmpPNEtime,PNEtmp.(OFFDATA.ChannelsFullName(i))(startPNE:endPNE),'-','LineWidth',1,'Marker','none','ShowBaselin','off')
    hold on
    
    
    %%%%Select which OFF periods to plot (duration/interval criteria)
    tmpOFFStarts=find(OFFDATA.StartOP(startPNE:endPNE,i));
    tmpOFFEnds=find(OFFDATA.EndOP(startPNE:endPNE,i));
    if ~isempty(tmpOFFStarts) 
        %Account for end section problem
        if tmpOFFStarts(1)>tmpOFFEnds(1)
            tmpOFFEnds=tmpOFFEnds(2:end);
        end
        if tmpOFFEnds(end)<tmpOFFStarts(end)
            tmpOFFStarts=tmpOFFStarts(1:end-1);
        end
        %Remove short gaps
        chooseON=tmpOFFStarts(2:end)-tmpOFFEnds(1:end-1)-1<(OFFDATA.PNEfs/1000*str2num(get(findobj('Tag','maxInterVal'),'String')));
        chooseONstart=~[0;chooseON];
        chooseONend=~[chooseON;0];
        tmpOFFStarts=tmpOFFStarts(chooseONstart);
        tmpOFFEnds=tmpOFFEnds(chooseONend);
        
        %Remove short off-periods
        chooseOFF=tmpOFFEnds-tmpOFFStarts+1>(OFFDATA.PNEfs/1000*str2num(get(findobj('Tag','minDurVal'),'String')));
        tmpOFFStarts=tmpOFFStarts(chooseOFF);
        tmpOFFEnds=tmpOFFEnds(chooseOFF);
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
    
    if i>1
        set(gca,'Xcolor','none')
    end
    
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

drawOFFP

end

function shiftDur(source,~)
set(findobj('Tag','minDurVal'),'String',num2str(source.Value))

drawOFFP
end

function shiftInter(source,~)
set(findobj('Tag','maxInterVal'),'String',num2str(source.Value))

drawOFFP
end

end
%clear PNEtmp
%clear LFPtmp
%%%% FAST LOADING EXAMPLE SCRIPT
%save('aa.mat','aa','-v7.3')
%exampleObject = matfile('D:\DPhil\Off-period detection\Matlab scripts\aa.mat');
%PNEtmp=exampleObject.aa(1,1:1000);