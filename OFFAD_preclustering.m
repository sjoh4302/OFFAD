function [OFFDATA,test]=OFFAD_preclustering(importDataVar)

set(findobj('Tag','OFFAD_IMPORT'),'Visible','Off')

%Plotting
g.Preclustering = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576/2 432], ...
	'PaperUnits','points', ...
	'name', 'OFFAD (OFF_period Automated Detection) - Pre-clustering', ... 
	'numbertitle', 'off', ...
	'Position',[200 100 800 450], ...
    'Toolbar','none',...
    'Menubar','none',...
	'Tag','OFFAD_PRECLUSTER',...
    'CloseRequestFcn',['set(findobj(''Tag'',''OFFAD_IMPORT''),''Visible'',''On'');'...
    'close(findobj(''Tag'',''OFFAD_PRECLUSTER''))']);
     
%Create empty directory to hold clustering information
OFFDATA=[];
%Add input variable from OFFAD_importdata
importDataVar=importDataVar(end:-1:1);
OFFDATA.datasetname=importDataVar(2).String;
OFFDATA.VSpathin=importDataVar(4).String;
OFFDATA.epochLen=double(string(importDataVar(7).String));
OFFDATA.PNEpathin=importDataVar(9).String;
OFFDATA.PNEfs=double(string(importDataVar(12).String));
OFFDATA.PNEunit=importDataVar(14).String(importDataVar(14).Value);
OFFDATA.LFPpathin=importDataVar(16).String;
OFFDATA.LFPfs=double(string(importDataVar(19).String));
OFFDATA.LFPunit=importDataVar(21).String(importDataVar(21).Value);
OFFDATA.FiltLFP=double(importDataVar(23).Value);
OFFDATA.clustEval=[importDataVar(25).String(importDataVar(25).Value)];
OFFDATA.ignoreChannels=double(string(regexp(importDataVar(27).String,'\d*','match')));
OFFDATA.clustVar1Select=double(string(importDataVar(29).Value)); %1 = amplitude, 2=power
OFFDATA.clustVar1Smooth=double(string(importDataVar(31).String));
OFFDATA.clustVar2Select=double(string(importDataVar(32).Value)); %1 = amplitude, 2=power
OFFDATA.clustVar2Smooth=double(string(importDataVar(34).String));
OFFDATA.percSamp=double(string(importDataVar(36).String)); %percentage of signal to sample for cluster thresholding
exampleObject = matfile(OFFDATA.PNEpathin);
channelNums = who(exampleObject);
channelNums = channelNums(find(contains(channelNums,'ch','IgnoreCase',true)));
%Sort channels
numOnly=double(string(regexp(string(channelNums),'\d*','match')));
[chanOrderNums,chanOrder]=sort(numOnly);
%Remove selected channels
if ~isempty(OFFDATA.ignoreChannels)
    chanOrder=(chanOrder(~abs(sum(chanOrderNums==OFFDATA.ignoreChannels,2))));
end
OFFDATA.ChannelsFullName=string(channelNums(chanOrder));
OFFDATA.Channels=double(string(cellfun(@(X) regexp(X,'\d*','match'),channelNums(chanOrder),'UniformOutput',false)));
OFFDATA.OptimalK=zeros(length(OFFDATA.Channels),1);
OFFDATA.BaselineAmp=nan(length(OFFDATA.Channels),1);
for nameChan=1:length(OFFDATA.ChannelsFullName)
    OFFDATA.GMModels.(OFFDATA.ChannelsFullName(nameChan))=[];
end
%% Initialise vigilance state data
% cutting off one 4 second epoch at the beginning and end of each episode to exclude transitional states
stage='nr';
%load vigilance state information for this animal and recording date ('nr' are all artefact free NREM epochs)
vigState = load(OFFDATA.VSpathin,'-mat',stage);
vigState = vigState.(stage);

% find NREM episodes
endEpi=find(diff(vigState)>1); %find last epoch of each episode
endEpi=[endEpi;length(vigState)]; %add final episode end
startEpi=find(diff(vigState)>1)+1; %find first epoch of each episode
startEpi=[1;startEpi]; %add first episode start
numepi=length(startEpi); %number of NREM episodes
global cleanepochs
global numepochs

cleanepochs=[];

%loop going through all NREM episodes 
for ep = 1:numepi
    
    Startepoch=vigState(startEpi(ep)); %find start epoch of this NREM episode
    Endepoch=vigState(endEpi(ep)); %find second last epoch of this NREM episode
    
    if Endepoch-Startepoch<2 %if NREM episode 2 epochs or shorter, there'll be no signal because we're cutting the first and last epoch (i.e. state transitions)
        continue
    else %if NREM episode at least 3 epochs, find bins corresponding to the start of the second epoch and end of the second last epoch
        cleanepochs=[cleanepochs; Startepoch Endepoch];
    end
end

numepochs=length(cleanepochs);


%%% Create button to return to import screen
uicontrol(g.Preclustering,'Style', 'pushbutton','String','Return',...
    'FontWeight','bold','FontSize',12,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Position',[0 0 0.08 0.1],'Callback','uiresume');

%%% Create button to continue with clustering
done=uicontrol(g.Preclustering,'Style', 'pushbutton','String','Done',...
    'FontWeight','bold','FontSize',12,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Tag','DoneButton',...
    'UserData',0,...
    'Position',[0.92 0 0.08 0.1],'Callback','set(findobj(''Tag'',''DoneButton''),''UserData'',1),uiresume');

%%% Create button to re-cluster channel
uicontrol(g.Preclustering,'Style', 'pushbutton','String','Recluster',...
    'FontWeight','bold','FontSize',12,...
    'BackgroundColor',[0.3 0.8 0.8],'Units','normalized',...
    'Position',[0.46 0 0.08 0.08],'Callback',@plotClustNum);

%%% Create list to scroll through channels
uicontrol(g.Preclustering,'Style', 'popupmenu','String',OFFDATA.ChannelsFullName,...
    'FontWeight','bold','FontSize',18,...
    'Units','normalized',...
    'Position',[0.2 0.95 0.6 0.05],...
    'CreateFcn',@plotClustNum,...
    'Callback',@plotClustNum);

 
uiwait
test=get(findobj('Tag','DoneButton'),'UserData');
set(findobj('Tag','OFFAD_IMPORT'),'Visible','On')
delete(findobj('Tag','OFFAD_PRECLUSTER'))


function plotClustNum(src,~)
       axesHandles=get(findobj('Tag','OFFAD_PRECLUSTER'),'Children');
       try
           cla(axesHandles(5));
           cla(axesHandles(6));
       end
       drawnow
 %%% load pNe signal
        channelString=axesHandles(1).String;
        channelValue=axesHandles(1).Value;
        if iscell(channelString)==0
            channelString={channelString};
        end
        PNE = load(OFFDATA.PNEpathin,'-mat',string(channelString(channelValue)));
        PNE = PNE.(string(channelString(channelValue)));
        
        if OFFDATA.PNEunit=="uV*1000000" %Convert to uV
            PNE=PNE/1000000;
        end
         
        PNE = int16(PNE);
        
        if length(unique(PNE))<4
            warning('Wrong amplitude units specified')
            return
        end
        
        PNE = abs(PNE); %take absolute values
        
       
                    
        %%% collect pNe signal for all episodes of chosen vigilance state
        vsPNE=NaN(1,length(PNE)); %NaN vectors to be filled with pNe signal
        vsPNEtime=NaN(1,length(PNE)); %NaN vectors to be filled with recording time values
         
        %loop going through all NREM epochs (except for first and last)
        for ep = 1:numepochs
            
            Startepoch=cleanepochs(ep,1); %start of respective epoch
            Endepoch=cleanepochs(ep,2); %end of respective epoch
            
            %fill NaN vectors with pNe signal for this NREM episode
            StartBin=ceil(Startepoch*OFFDATA.epochLen*OFFDATA.PNEfs);
            EndBin=floor((Endepoch-1)*OFFDATA.epochLen*OFFDATA.PNEfs);
            try
                vsPNE(StartBin:EndBin)=PNE(StartBin:EndBin);
                vsPNEtime(StartBin:EndBin)=StartBin:EndBin;
            end
            
        end
        %concatenate the signal from all selected NREM episodes to get rid of gaps
        vsPNE=vsPNE(~isnan(vsPNE));
        vsPNEtime=vsPNEtime(~isnan(vsPNEtime));
        clear PNE
        
%end

%% Data processing
%%%Cluster 1
if OFFDATA.clustVar1Select==1
    winSize=(OFFDATA.clustVar1Smooth*2)+1;
    clusterVar1=conv(vsPNE,gausswin(winSize));
    if mod(winSize,2)==0
        clusterVar1=clusterVar1(winSize/2:end-winSize/2);
    else
        clusterVar1=clusterVar1(ceil(winSize/2):end-floor(winSize/2));
    end
    clear winSize
elseif OFFDATA.clustVar1Select==2
    percOverlap=10; %overlap between segments
    LB_freq=40; %find signal power from LB to fs/2
    [clusterVar1,F] = offPeriodScalogram(vsPNE,OFFDATA.PNEfs,percOverlap,LB_freq);
    clear LB_freq percOverlap F
end

%%%Cluster 2
if OFFDATA.clustVar2Select==1
    winSize=(OFFDATA.clustVar2Smooth*2)+1;
    clusterVar2=conv(vsPNE,gausswin(winSize));
    if mod(winSize,2)==0
        clusterVar2=clusterVar2(winSize/2:end-winSize/2);
    else
        clusterVar2=clusterVar2(ceil(winSize/2):end-floor(winSize/2));
    end
    clear winSize
elseif OFFDATA.clustVar2Select==2
    percOverlap=10; %overlap between segments
    LB_freq=40; %find signal power from LB to fs/2
    [clusterVar2,F] = offPeriodScalogram(vsPNE,OFFDATA.PNEfs,percOverlap,LB_freq);
    clear LB_freq percOverlap F
end

%% Cluster sample data 


    randPoints=randi(length(clusterVar1),round(length(vigState)*4*OFFDATA.PNEfs*(1/1000)),1);
    sampCluster=[clusterVar1(randPoints)',clusterVar2(randPoints)'];
    allIDX=[];
    %allIDX2=[];
    %%%DONT SELECT SHARED COVARIANCE (fails to find off period cluster)
    for i = 1:8
        options = statset('MaxIter',1000,'TolFun',1e-5);
        try
            GMModels = fitgmdist(sampCluster,i,'Replicates',5,'Options',options);
            allIDX(:,i)=cluster(GMModels,sampCluster);
        catch
            warning('Failed to cluster channel')
            return %If clustering fails, cancel
        end
    end
   % if size(allIDX,2)==8 %Only proceed if clustering succesful
    eva = evalclusters(sampCluster,allIDX,string(OFFDATA.clustEval));
    OFFDATA.OptimalK(channelValue)=eva.OptimalK;
    GMModel = fitgmdist(sampCluster,eva.OptimalK,'Replicates',5,'Options',options);
   
%% Plot cluster examples
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
    subplot(1,2,2)
    d = 300; % Grid density
    x1 = linspace(min(clusterVar1), prctile(clusterVar1,99), d);
    x2 = linspace(min(clusterVar2), prctile(clusterVar2,99), d);
    [x1grid,x2grid] = meshgrid(x1,x2);
    X0 = [x1grid(:) x2grid(:)];
    mahalDist = mahal(GMModel,X0);
    threshold = sqrt(chi2inv(0.99,2));
    hold on
    [~,clustOrder]=sort(GMModel.mu(:,1)); 
    for m = 1:size(mahalDist,2)
        currentClust=clustOrder(m);
        idx = mahalDist(:,currentClust)<=threshold;
        if m==1
            h2 = plot(X0(idx,1),X0(idx,2),'r.','MarkerSize',1);
        else
           h2 = plot(X0(idx,1),X0(idx,2),'b.','MarkerSize',1);
        end
        uistack(h2,'bottom');
    end 
    xlabel('Clustering variable 1')
    ylabel('Clustering variable 2')
    xlmts=get(gca,'XLim');
    ylmts=get(gca,'YLim');
    text(xlmts(1)+0.05*(xlmts(2)-xlmts(1)),ylmts(1)+0.9*(ylmts(2)-ylmts(1)),...
    {['Clust Eval: ',char(string(OFFDATA.clustEval))],['Optimal K: ',char(string(eva.OptimalK))]},...
    'FontSize',14)
    
     
    clusterAxes=gca;
    subplot(1,2,1)
    Var1Min=clusterAxes.XLim(1);
    Var1Max=clusterAxes.XLim(2);
    Var2Min=clusterAxes.YLim(1);
    Var2Max=clusterAxes.YLim(2);
    nhistBins=200;
    Var1edges=Var1Min:(Var1Max-Var1Min)/nhistBins:Var1Max;
    Var2edges=Var2Min:(Var2Max-Var2Min)/nhistBins:Var2Max;
    [N,Var1edges,Var2edges] = histcounts2(clusterVar1,clusterVar2,Var1edges,Var2edges);
    N(N==0)= max(max(N));
    imagesc('XData',[Var1Min,Var1Max],'YData',[Var2Min,Var2Max],'CData',N')
    xlim([Var1Min,Var1Max])
    ylim([Var2Min,Var2Max])
    c = hot;
    %c = flipud(c);
    colormap(c);
    xlabel('Clustering variable 1')
    ylabel('Clustering variable 2')
   % caxis([0 2000])
    %set(gca,'ColorScale','log')
    %set(gca,'CLim',[0.1 2000])
    
%     figure
%     IDX = cluster(GMModel,sampCluster);
%     gscatter(sampCluster(:,1),sampCluster(:,2),IDX)
%     xlim([clusterAxes.XLim(1),clusterAxes.XLim(2)])
%     ylim([clusterAxes.YLim(1),clusterAxes.YLim(2)])


end

end
