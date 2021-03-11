function [OFFDATA]=OFFAD_scrollSave(OFFDATA,minDur,maxInt);
g.Save = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576/2 432], ...
	'PaperUnits','points', ...
	'name', 'OFFAD (OFF_period Automated Detection) - Clustering', ... 
	'numbertitle', 'off', ...
	'Position',[200 300 700 100], ...
    'Toolbar','none',...
    'Menubar','none',...
	'Tag','OFFAD_SCROLLSAVE');
uicontrol(g.Save,'Style', 'text','String','Warning: Saving adjusted OFF periods',...
    'FontWeight','bold','FontSize',22,...
    'Units','normalized',...
    'Position',[0.2 0.2 0.6 0.6]);

drawnow 

OFFDATA.StartOPadjusted=sparse(repmat(logical(0),length(OFFDATA.StartOP),length(OFFDATA.ChannelsFullName)));
OFFDATA.EndOPadjusted=sparse(repmat(logical(0),length(OFFDATA.StartOP),length(OFFDATA.ChannelsFullName)));
OFFDATA.AllOPadjusted=sparse(repmat(logical(0),length(OFFDATA.StartOP),length(OFFDATA.ChannelsFullName)));

numChan=length(OFFDATA.ChannelsFullName);
for i = 1:numChan
        adjOFFStarts=find(OFFDATA.StartOP(:,i));
        adjOFFEnds=find(OFFDATA.EndOP(:,i));
        %Remove short gaps
        chooseON=adjOFFStarts(2:end)-adjOFFEnds(1:end-1)-1<(OFFDATA.PNEfs/1000*str2num(maxInt));
        chooseONstart=~[0;chooseON];
        chooseONend=~[chooseON;0];
        adjOFFStarts=adjOFFStarts(chooseONstart);
        adjOFFEnds=adjOFFEnds(chooseONend);
        
        %Remove short off-periods
        chooseOFF=adjOFFEnds-adjOFFStarts+1>(OFFDATA.PNEfs/1000*str2num(minDur));
        adjOFFStarts=adjOFFStarts(chooseOFF);
        adjOFFEnds=adjOFFEnds(chooseOFF);
       
        
        adjOFFall=[];
        for k = 1:length(adjOFFStarts)
            adjOFFall=[adjOFFall,[adjOFFStarts(k):adjOFFEnds(k)]];
        end
        
        %Store OFF duration threshold
        OFFDATA.OFFthresh=str2num(minDur);
        
        %Store ON duration threshold
        OFFDATA.ONthresh=str2num(maxInt);
        
        %Store START OFF-P data
        OFFDATA.StartOPadjusted(:,numChan)=sparse(adjOFFStarts,1,logical(1),length(OFFDATA.StartOPadjusted),1,length(adjOFFStarts));

        %Store END OFF-P data
        OFFDATA.EndOPadjusted(:,numChan)=sparse(adjOFFEnds,1,logical(1),length(OFFDATA.StartOPadjusted),1,length(adjOFFEnds));

        %Store ALL OFF-P data
        OFFDATA.AllOPadjusted(:,numChan)=sparse(adjOFFall,1,logical(1),length(OFFDATA.StartOPadjusted),1,length(adjOFFall));

end        
        
close(findobj('Tag','OFFAD_SCROLLSAVE'))

end