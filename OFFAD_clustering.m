function [OFFDATA]=OFFAD_clustering(OFFDATA)
functionTimer = tic;

%Plotting
g.Clustering = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576/2 432], ...
	'PaperUnits','points', ...
	'name', 'OFFAD (OFF_period Automated Detection) - Clustering', ... 
	'numbertitle', 'off', ...
	'Position',[200 300 700 100], ...
    'Toolbar','none',...
    'Menubar','none',...
	'Tag','OFFAD_CLUSTER');
uicontrol(g.Clustering,'Style', 'text','String','Warning: Clustering in progress',...
    'FontWeight','bold','FontSize',22,...
    'Units','normalized',...
    'Position',[0.2 0.2 0.6 0.6]);

drawnow 


%%%%% Perform clustering
%% collect bin number for start and end of each NREM epoch
% cutting off one 4 second epoch at the beginning and end of each episode to exclude transitional states

%load vigilance state information for this animal and recording date ('nr' are all artefact free NREM epochs)
vigState = load(OFFDATA.VSpathin,'-mat','nr');
vigState = vigState.('nr');

% find NREM episodes
endEpi=find(diff(vigState)>1); %find last epoch of each episode
endEpi=[endEpi;length(vigState)]; %add final episode end
startEpi=find(diff(vigState)>1)+1; %find first epoch of each episode
startEpi=[1;startEpi]; %add first episode start
numepi=length(startEpi); %number of NREM episodes
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


%% Extract raw pNe data

%{
1) Convert it into absolute values
2) Extract and concatenate NREM episodes only
%}

exampleObject = matfile(OFFDATA.PNEpathin);

PNElength=size(exampleObject,OFFDATA.ChannelsFullName(1),2);
OFFDATA.StartOP=sparse(repmat(logical(0),PNElength,length(OFFDATA.ChannelsFullName)));
OFFDATA.EndOP=sparse(repmat(logical(0),PNElength,length(OFFDATA.ChannelsFullName)));
OFFDATA.AllOP=sparse(repmat(logical(0),PNElength,length(OFFDATA.ChannelsFullName)));
clear exampleObject 

for chanNum = 1:length(OFFDATA.ChannelsFullName)
    %Start channel timer
    channelTimer=tic;
    
    chan=OFFDATA.ChannelsFullName(chanNum);    
    display(['Clustering ',char(chan)])

    %skip channel if it's not a good channel (was excluded on visual inspection due to bad signals)
            
        %%% load pNe signal
        PNE = load(OFFDATA.PNEpathin,'-mat',chan);
        PNE = PNE.(chan);
        PNE = int16(PNE);
        
        PNE = abs(PNE); %take absolute values 
        
        if OFFDATA.PNEunit=="V" %Convert to uV
            PNE=PNE/1000000;
        end
                        
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
        
        
    
%end



%% Data processing
%%%Cluster 1
if OFFDATA.clustVar1Select==1
    winSize=(OFFDATA.clustVar1Smooth*2)+1;
    clusterVar1=conv(vsPNE,gausswin(winSize));
    clusterVar1_full=conv(PNE,gausswin(winSize));
    if mod(winSize,2)==0
        clusterVar1=clusterVar1(winSize/2:end-winSize/2);
        clusterVar1_full=clusterVar1_full(winSize/2:end-winSize/2);
    else
        clusterVar1=clusterVar1(ceil(winSize/2):end-floor(winSize/2));
        clusterVar1_full=clusterVar1_full(ceil(winSize/2):end-floor(winSize/2));
    end
    clear winSize
elseif OFFDATA.clustVar1Select==2
    percOverlap=10; %overlap between segments
    LB_freq=40; %find signal power from LB to fs/2
    [clusterVar1,F] = offPeriodScalogram(vsPNE,OFFDATA.PNEfs,percOverlap,LB_freq);
    [clusterVar1_full,F] = offPeriodScalogram(PNE,OFFDATA.PNEfs,percOverlap,LB_freq);
    clear LB_freq percOverlap F
end

%%%Cluster 2
if OFFDATA.clustVar2Select==1
    winSize=(OFFDATA.clustVar2Smooth*2)+1;
    clusterVar2=conv(vsPNE,gausswin(winSize));
     clusterVar2_full=conv(PNE,gausswin(winSize));
    if mod(winSize,2)==0
        clusterVar2=clusterVar2(winSize/2:end-winSize/2);
        clusterVar2_full=clusterVar2_full(winSize/2:end-winSize/2);
    else
        clusterVar2=clusterVar2(ceil(winSize/2):end-floor(winSize/2));
        clusterVar2_full=clusterVar2_full(ceil(winSize/2):end-floor(winSize/2));
    end
    clear winSize
elseif OFFDATA.clustVar2Select==2
    percOverlap=10; %overlap between segments
    LB_freq=40; %find signal power from LB to fs/2
    [clusterVar2,F] = offPeriodScalogram(vsPNE,OFFDATA.PNEfs,percOverlap,LB_freq);
    [clusterVar2_full,F] = offPeriodScalogram(PNE,OFFDATA.PNEfs,percOverlap,LB_freq);
    clear LB_freq percOverlap F
end

clear PNE
clear vsPNE

%% Gaussian clustering

if OFFDATA.OptimalK(chanNum)==0
    randPoints=randi(length(clusterVar1),10000,1);
    sampCluster=[clusterVar1(randPoints)',clusterVar2(randPoints)'];
    allIDX=[];
    %allIDX2=[];
    %%%DONT SELECT SHARED COVARIANCE (fails to find off period cluster)
    for i = 1:8
        options = statset('MaxIter',1000,'TolFun',1e-5);
        GMModels = fitgmdist(sampCluster,i,'Replicates',5,'Options',options);
        allIDX(:,i)=cluster(GMModels,sampCluster);
        %allIDX2(:,i)=cluster(GMModels,clusterData2);
    end
    eva = evalclusters(sampCluster,allIDX,string(OFFDATA.clustEval));
    OFFDATA.OptimalK(chanNum)=eva.OptimalK;
    clear allIDX sampCluster randPoints
end

allData=[clusterVar1_full',clusterVar2_full'];
randPoints=randi(length(clusterVar1),round(length(clusterVar1)*(OFFDATA.percSamp/100)),1);
sampData=[clusterVar1(randPoints)',clusterVar2(randPoints)'];
options = statset('MaxIter',1000,'TolFun',1e-5);
GMModel = fitgmdist(sampData,OFFDATA.OptimalK(chanNum),'Replicates',20,'Options',options);
IDX = cluster(GMModel,allData);


%% Find off periods 
%profile on
allPoints=[1:1:PNElength];

[~,offCluster]=min(GMModel.mu(:,1));
OFF_clust_points=allPoints(IDX==offCluster);

%%%%% Find all OFF periods
OFFgaps=find(diff(OFF_clust_points)>1); %find last epoch of each episode
numOFF=length(OFFgaps); %number of OFF periods
%loop going through all OFF periods
for ep = 1:numOFF+1
    if ep==1
        StartOFF(ep)=OFF_clust_points(1); %find start point of this OFF period
        EndOFF(ep)=OFF_clust_points(OFFgaps(ep)); %find last point of this OFF period
    elseif ep==numOFF+1
        StartOFF(ep)=OFF_clust_points(OFFgaps(ep-1)+1); %find start point of this OFF period
        EndOFF(ep)=OFF_clust_points(end); %find last point of this OFF period
    else
        StartOFF(ep)=OFF_clust_points(OFFgaps(ep-1)+1); %find start point of this OFF period
        EndOFF(ep)=OFF_clust_points(OFFgaps(ep)); %find last point of this OFF period
    end
    
   
end

%Store START OFF-P data
OFFDATA.StartOP(:,chanNum)=sparse(StartOFF,1,logical(1),length(OFFDATA.StartOP),1,length(StartOFF));

%Store END OFF-P data
OFFDATA.EndOP(:,chanNum)=sparse(EndOFF,1,logical(1),length(OFFDATA.StartOP),1,length(EndOFF));

%Store ALL OFF-P data
OFFDATA.AllOP(:,chanNum)=sparse(OFF_clust_points,1,logical(1),length(OFFDATA.StartOP),1,length(OFF_clust_points));


clear OFFperiod OFF_clust_points ON_clust_points clusterVar1 clusterVar2 ...
      clusterIDX_total clustVar1ThreshScaled clustVar2ThreshScaled StartOFF EndOFF

%Display total channel clustering time
channelTime=toc(channelTimer);
display(['Channel clustering time = ',char(string(channelTime)),' sec'])
end

%Display total function clustering time
functionTime=toc(functionTimer);
display(['Total clustering time = ',char(string(functionTime)),' sec'])

%Return to main menu    
close(findobj('Tag','OFFAD_IMPORT'));
OFF_AD('redraw');
set(findobj('Tag','OFFAD_CLUSTER'),'Visible','Off');


%n=10;
%S = spalloc(n,n,3*n);
%for j = 1:n
%    ind = [max(j-1,1) j min(j+1,n)]
%    S(:,j) = sparse(ind,1,round(rand(3,1)),n,1,3);
%end