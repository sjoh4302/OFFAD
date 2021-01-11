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
LFPtmp=load(OFFDATA.LFPpathin,OFFDATA.ChannelsFullName{:});

PNEexampleObject = matfile(OFFDATA.PNEpathin);
PNElength=size(PNEexampleObject,OFFDATA.ChannelsFullName(1),2);

LFPexampleObject = matfile(OFFDATA.LFPpathin);
LFPlength=size(LFPexampleObject,OFFDATA.ChannelsFullName(1),2);


uicontrol(g.Scroll,'Style', 'text','String','Seconds',...
    'FontWeight','bold','FontSize',10,...
    'Units','normalized',...
    'Position',[0.54 0.01 0.08 0.04]);

uicontrol(g.Scroll,'Style', 'edit','String','0',...
    'FontWeight','bold','FontSize',10,...
    'Units','normalized',...
    'Position',[0.46 0.01 0.08 0.04],'CreateFcn',@drawOFFP,...
    'Callback',@drawOFFP);
 

function drawOFFP(source,~)
startSec=str2num(source.String);
startPNE=round(startSec/(1/OFFDATA.PNEfs))+1;
endPNE=round(startSec/(1/OFFDATA.PNEfs))+2000;

tmpPNEtime=[startPNE/OFFDATA.PNEfs:1/OFFDATA.PNEfs:endPNE/OFFDATA.PNEfs];
startLFP=round(mod(startPNE/OFFDATA.PNEfs,1/OFFDATA.LFPfs)+startPNE/OFFDATA.PNEfs/(1/OFFDATA.LFPfs));
endLFP=round(mod(endPNE/OFFDATA.PNEfs,1/OFFDATA.LFPfs)+endPNE/OFFDATA.PNEfs/(1/OFFDATA.LFPfs));
tmpLFPtime=[startLFP/OFFDATA.LFPfs:1/OFFDATA.LFPfs:endLFP/OFFDATA.LFPfs];


numChan=length(OFFDATA.Channels);
subplot2=linspace(0.1,0.95,numChan+1);
subplot4=(subplot2(2)-subplot2(1))-(subplot2(2)-subplot2(1))/10;
for i = 1:numChan
    %PNEtmp=load(OFFDATA.PNEpathin,OFFDATA.ChannelsFullName(i));
    %PNEtmp=PNEtmp.(OFFDATA.ChannelsFullName(i));
    subplot('Position',[0.06 subplot2(i) 0.88 subplot4]);
    yyaxis left
    stem(tmpPNEtime,PNEtmp.(OFFDATA.ChannelsFullName(i))(startPNE:endPNE),'-','LineWidth',1,'Marker','none','ShowBaselin','off')
    hold on
    stem(tmpPNEtime(find(OFFDATA.nr.AllOP(startPNE:endPNE,i))),PNEtmp.(OFFDATA.ChannelsFullName(i))(startPNE-1+find(OFFDATA.nr.AllOP(startPNE:endPNE,i))),'-','LineWidth',1,'Marker','none','ShowBaselin','off','Color','red','LineStyle','-')
    ylim([-150 150])
    
    %LFPtmp=load(OFFDATA.LFPpathin,OFFDATA.ChannelsFullName(i));
    %LFPtmp=LFPtmp.(OFFDATA.ChannelsFullName(i));
    yyaxis right
    plot(tmpLFPtime,LFPtmp.(OFFDATA.ChannelsFullName(i))(startLFP:endLFP),'-')
    xlim([tmpPNEtime(1) tmpPNEtime(end)])
    set(gca,'Box','off')
    if i>1
        set(gca,'Xcolor','none')
    end
    
end

end

end
%clear PNEtmp
%clear LFPtmp
%%%% FAST LOADING EXAMPLE SCRIPT
%save('aa.mat','aa','-v7.3')
%exampleObject = matfile('D:\DPhil\Off-period detection\Matlab scripts\aa.mat');
%PNEtmp=exampleObject.aa(1,1:1000);