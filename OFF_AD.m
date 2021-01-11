function OFF_AD( arg )
if nargin < 1
    hh = findobj('tag', 'OFFAD');
    if ~isempty(hh)
        disp('OFFAD warning: there can be only one OFFAD window, closing old one');
        close(hh);  
    end
    
    evalin('base', 'global OFFDATA');
    global OFFDATA
     OFFDATA=[];
    
g.Main = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'name', 'OFFAD (OFF_period Automated Detection)', ... 
	'numbertitle', 'off', ...
	'Position',[200 100 500 400], ...
	'Toolbar','none',...
    'Menubar','none',...
    'Tag','OFFAD');

uicontrol(g.Main,'Style', 'pushbutton','String','Open new study',...
    'FontWeight','bold','FontSize',20,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Position',[0.2 0.55 0.6 0.35],...
    'Tag','M_IMPORT',...
    'Callback',['[OFFDATA] = OFFAD_importdata;'...
    'waitfor(findobj(''Tag'',''OFFAD''),''Visible'');' 'OFF_AD(''redraw'');']);

uicontrol(g.Main,'Style', 'pushbutton','String','Load  previous study',...
    'FontWeight','bold','FontSize',20,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Tag','M_LOAD',...
    'Position',[0.2 0.1 0.6 0.35],'Callback',['uiload;' 'OFF_AD(''redraw'');']);

else
    
    set(findobj('Tag','M_IMPORT'),'Visible','Off')
    set(findobj('Tag','M_LOAD'),'Visible','Off')
    set(findobj('Tag','OFFAD'),'Visible','On')
    
    g.Main=findobj('Tag','OFFAD');
    
    global OFFDATA
    
    uicontrol(g.Main,'Style', 'text','String',['Dataset: ' OFFDATA.datasetname],...
    'FontWeight','bold','FontSize',18,...
    'Units','normalized',...
    'Position',[0.2 0.85 0.6 0.1]);
    
    uicontrol(g.Main,'Style', 'pushbutton','String','Scroll through off-periods',...
    'FontWeight','bold','FontSize',18,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Position',[0.2 0.6 0.6 0.15],...
    'Tag','M_VIEWER',...
    'Callback','OFFAD_scroll(OFFDATA);');

    uicontrol(g.Main,'Style', 'pushbutton','String','View channel statistics',...
    'FontWeight','bold','FontSize',18,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Position',[0.2 0.4 0.6 0.15],...
    'Tag','M_CHANSTAT',...
    'Callback','OFFAD_channelstats(OFFDATA);');

    uicontrol(g.Main,'Style', 'pushbutton','String','Save to current folder',...
    'FontWeight','bold','FontSize',18,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Position',[0.2 0.2 0.6 0.15],...
    'Tag','M_SAVECURR',...
    'Callback','save([OFFDATA.datasetname,''.mat''],''OFFDATA'')');

    uicontrol(g.Main,'Style', 'pushbutton','String','Save to new folder',...
    'FontWeight','bold','FontSize',18,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Position',[0.2 0.02 0.6 0.15],...
    'Tag','M_SAVECURR',...
    'Callback','uisave(''OFFDATA'',OFFDATA.datasetname)');
 
 
end

 % global OFFDATA_var

%global OFFDATA
%global OFFDATA_var

%if isempty(OFFDATA), OFFDATA = []; end
%if isempty(OFFDATA_var), OFFDATA_var = []; end

