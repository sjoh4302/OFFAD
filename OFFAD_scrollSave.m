function [OFFDATA]=OFFAD_scrollSave(OFFDATA,minDur,maxInt);
%
% Scroll saving page: Save adjusted OFF periods with interval and duration thresholds applied
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
    
% Create figure for saving window
g.Save = figure('name', 'OFFAD (OFF_period Automated Detection) - Save thresholds', ... 
	'numbertitle', 'off', ...
	'Position',[screensize(3)/6,3*screensize(4)/8,4*screensize(3)/6,screensize(4)/4], ...
    'Toolbar','none',...
    'Menubar','none',...
	'Tag','OFFAD_SCROLLSAVE');
uicontrol(g.Save,'Style', 'text','String','Warning: Storing adjusted OFF periods',...
    'FontWeight','bold','FontSize',22,...
    'Units','normalized',...
    'Position',[screensize(3)/6,3*screensize(4)/8,4*screensize(3)/6,screensize(4)/4])
drawnow

% Initalise/reset adjusted OFF period START,END and ALL data point matrices
OFFDATA.StartOPadjusted=sparse(repmat(logical(0),length(OFFDATA.StartOP),length(OFFDATA.ChannelsFullName)));
OFFDATA.EndOPadjusted=sparse(repmat(logical(0),length(OFFDATA.StartOP),length(OFFDATA.ChannelsFullName)));
OFFDATA.AllOPadjusted=sparse(repmat(logical(0),length(OFFDATA.StartOP),length(OFFDATA.ChannelsFullName)));

numChan=length(OFFDATA.ChannelsFullName);
for i = 1:numChan

    % Find new OFF period start and end times 
    adjOFFStarts=find(OFFDATA.StartOP(:,i));
    adjOFFEnds=find(OFFDATA.EndOP(:,i));

    % Check that channel contains OFF periods before applying thresholds
    if ~isempty(adjOFFStarts)

        % Apply maximum tolerated intra-OFF period interval threshold
        chooseON=adjOFFStarts(2:end)-adjOFFEnds(1:end-1)-1<(OFFDATA.MUAfs/1000*maxInt(i));
        chooseONstart=~[0;chooseON];
        chooseONend=~[chooseON;0];
        adjOFFStarts=adjOFFStarts(chooseONstart);
        adjOFFEnds=adjOFFEnds(chooseONend);

        % Apply minimum OFF period duration threshold
        chooseOFF=adjOFFEnds-adjOFFStarts+1>(OFFDATA.MUAfs/1000*minDur(i));
        adjOFFStarts=adjOFFStarts(chooseOFF);
        adjOFFEnds=adjOFFEnds(chooseOFF);
    end

    adjOFFall=[];
    for k = 1:length(adjOFFStarts)
        adjOFFall=[adjOFFall,[adjOFFStarts(k):adjOFFEnds(k)]];
    end

    %Store OFF duration threshold
    OFFDATA.OFFthresh=minDur;

    %Store ON duration threshold
    OFFDATA.ONthresh=maxInt;

    %Store START OFF-P data
    OFFDATA.StartOPadjusted(:,i)=sparse(adjOFFStarts,1,logical(1),length(OFFDATA.StartOPadjusted),1,length(adjOFFStarts));

    %Store END OFF-P data
    OFFDATA.EndOPadjusted(:,i)=sparse(adjOFFEnds,1,logical(1),length(OFFDATA.StartOPadjusted),1,length(adjOFFEnds));

    %Store ALL OFF-P data
    OFFDATA.AllOPadjusted(:,i)=sparse(adjOFFall,1,logical(1),length(OFFDATA.StartOPadjusted),1,length(adjOFFall));

    clear adjOFFStarts adjOFFEnds chooseON chooseONstart chooseONend chooseOFF 
end        
        
close(findobj('Tag','OFFAD_SCROLLSAVE'))

end