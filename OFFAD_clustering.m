function [OFFDATA,OFFDATA_var]=OFFAD_clustering(importDataVar)

g.Import = findobj('tag', 'OFFAD_IMPORT');
set(g.Import,'Visible','off')

%Plotting
g.Clustering = figure('Units','points', ...
... %	'Colormap','gray', ...
	'PaperPosition',[18 180 576/2 432], ...
	'PaperUnits','points', ...
	'name', 'OFFAD (OFF_period Automated Detection)', ... 
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

% Extract import information from OFFAD_importdata
OFFDATA_var=[];
importDataVar=importDataVar(end:-1:1);
OFFDATA_var.datasetname=importDataVar(2).String;
OFFDATA_var.VSpathin=importDataVar(4).String;
OFFDATA_var.epochLen=double(string(importDataVar(7).String));
OFFDATA_var.PNEpathin=importDataVar(9).String;
OFFDATA_var.PNEfs=double(string(importDataVar(12).String));
OFFDATA_var.PNEunits=importDataVar(14).String;
OFFDATA_var.LFPpathin=importDataVar(16).String;
OFFDATA_var.LFPfs=double(string(importDataVar(19).String));
OFFDATA_var.LFPunit=importDataVar(21).String;
OFFDATA_var.selectStages=[importDataVar(24).Value,...      
              importDataVar(26).Value,...
              importDataVar(28).Value];
OFFDATA_var.ignoreChannels=double(string(regexp(importDataVar(30).String,'\d*','match')));
OFFDATA_var.clustVar1Select=double(string(importDataVar(32).Value)); %1 = amplitude, 2=power
OFFDATA_var.clustVar1Smooth=double(string(importDataVar(34).String));
OFFDATA_var.clustVar2Select=double(string(importDataVar(35).Value)); %1 = amplitude, 2=power
OFFDATA_var.clustVar2Smooth=double(string(importDataVar(37).String));
OFFDATA_var.percSamp=double(string(importDataVar(39).String)); %percentage of signal to sample for cluster thresholding

%Create empty directory to hold clustering information
OFFDATA=[];

%%%%% Perform clustering
%Choose stage to analyse
allstages=["nr","r","w"];
stages=allstages(logical(OFFDATA_var.selectStages));
for stageNum=1:length(stages)
stage=stages(stageNum);
%% collect bin number for start and end of each NREM epoch
% cutting off one 4 second epoch at the beginning and end of each episode to exclude transitional states

%load vigilance state information for this animal and recording date ('nr' are all artefact free NREM epochs)
vigState = load(OFFDATA_var.VSpathin,'-mat',stage);
vigState = vigState.(stage);

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

exampleObject = matfile(OFFDATA_var.PNEpathin);
channelNums = who(exampleObject);
%Sort channels
numOnly=double(string(regexp(string(channelNums),'\d*','match')));
[chanOrderNums,chanOrder]=sort(numOnly);
%Remove selected channels
if ~isempty(OFFDATA_var.ignoreChannels)
    chanOrder=(chanOrder(~abs(sum(chanOrderNums==OFFDATA_var.ignoreChannels,2))));
end
clear exampleObject


for chanNum = 1:length(chanOrder)
    chan=string(channelNums(chanOrder(chanNum)));    
    
    %skip channel if it's not a good channel (was excluded on visual inspection due to bad signals)
            
        %%% load pNe signal
        PNE = load(OFFDATA_var.PNEpathin,'-mat',chan);
        PNE = PNE.(chan);
        
        PNE = abs(PNE); %take absolute values
        %PNE = PNE(1:1000000);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%REMOV!!!!!!!!!!!!!!!!!
        timePNE=[1/OFFDATA_var.PNEfs:1/OFFDATA_var.PNEfs:length(PNE)/OFFDATA_var.PNEfs];
        
        
        %%% collect pNe signal for all episodes of chosen vigilance state
        vsPNE=NaN(1,length(PNE)); %NaN vectors to be filled with pNe signal
        vsPNEtime=NaN(1,length(PNE)); %NaN vectors to be filled with recording time values
         
        %loop going through all NREM epochs (except for first and last)
        for ep = 1:numepochs
            
            Startepoch=cleanepochs(ep,1); %start of respective epoch
            Endepoch=cleanepochs(ep,2); %end of respective epoch
            
            %fill NaN vectors with pNe signal for this NREM episode
            StartBin=ceil(Startepoch*OFFDATA_var.epochLen*OFFDATA_var.PNEfs);
            EndBin=floor((Endepoch-1)*OFFDATA_var.epochLen*OFFDATA_var.PNEfs);
            try
                vsPNE(StartBin:EndBin)=PNE(StartBin:EndBin);
                vsPNEtime(StartBin:EndBin)=timePNE(StartBin:EndBin);
            end
            
        end
        %concatenate the signal from all selected NREM episodes to get rid of gaps
        vsPNE=vsPNE(~isnan(vsPNE));
        vsPNEtime=vsPNEtime(~isnan(vsPNEtime));
        clear PNE
        
    
%end



%% Data processing
%%%Cluster 1
if OFFDATA_var.clustVar1Select==1
    winSize=OFFDATA_var.clustVar1Smooth;
    clusterVar1=conv(vsPNE,gausswin(winSize));
    if mod(winSize,2)==0
        clusterVar1=clusterVar1(winSize/2:end-winSize/2);
    else
        clusterVar1=clusterVar1(ceil(winSize/2):floor(winSize/2));
    end
    clear winSize
elseif OFFDATA_var.clustVar1Select==2
    percOverlap=10; %overlap between segments
    LB_freq=40; %find signal power from LB to fs/2
    [clusterVar1,F] = offPeriodScalogram(vsPNE,OFFDATA_var.PNEfs,percOverlap,LB_freq);
    clear LB_freq percOverlap F
end

%%%Cluster 2
if OFFDATA_var.clustVar2Select==1
    winSize=OFFDATA_var.clustVar2Smooth;
    clusterVar2=conv(vsPNE,gausswin(winSize));
    if mod(winSize,2)==0
        clusterVar2=clusterVar2(winSize/2:end-winSize/2);
    else
        clusterVar2=clusterVar2(ceil(winSize/2):floor(winSize/2));
    end
    clear winSize
elseif OFFDATA_var.clustVar2Select==2
    percOverlap=10; %overlap between segments
    LB_freq=40; %find signal power from LB to fs/2
    [clusterVar2,F] = offPeriodScalogram(vsPNE,OFFDATA_var.PNEfs,percOverlap,LB_freq);
    clear LB_freq percOverlap F
end

%% Spectral clustering (3 clusters)
%Extracts 30s samples of recording and perfroms spectral clustering
sampleLength=30; %in seconds
totalSamp=floor(((length(vsPNE)/OFFDATA_var.PNEfs)/sampleLength));
sampSize=round(OFFDATA_var.PNEfs*sampleLength);
sampStarts=([0:1:totalSamp-1]*sampSize)+1;
numSamp=round(totalSamp*(OFFDATA_var.percSamp/100));
color = lines(3);
sampStarts=sampStarts(randsample(length(sampStarts),numSamp));

%Declare global variables for thresholding function
global X
global idx
global OFFclusterINDEX

for n = 1:numSamp 
startPoint=sampStarts(n);
display([char(chan),' ',char(stage),' Sample:',num2str(n),'/',num2str(numSamp)])

X=[clusterVar1(startPoint:startPoint+sampSize-1)',clusterVar2(startPoint:startPoint+sampSize-1)'];
Xscaled=[rescale(X(:,1)),rescale(X(:,2))];

idx = spectralcluster(Xscaled,3);

%%% Find thresholds from clusterings
[clustVar1ThreshScaled,OFFclusterINDEX]=min([max(Xscaled(idx==1,1)),max(Xscaled(idx==2,1)),max(Xscaled(idx==3,1))]);
clustVar1Thresh=X(find(Xscaled(:,1)==clustVar1ThreshScaled),1);
[clustVar2ThreshScaled]=min([max(Xscaled(idx==1,2)),max(Xscaled(idx==2,2)),max(Xscaled(idx==3,2))]);
clustVar2Thresh=X(find(Xscaled(:,2)==clustVar2ThreshScaled),2);

% Use fminsearch to find thresholds which minimse the loss function (false
% positive plus false negative
x0 = [clustVar1Thresh,clustVar2Thresh];
[sampleThresh,fval] = fminsearch(@clusterThresholds,x0);
allSampleThresholds(n,:)=sampleThresh;
end

%% Find off periods using clustered thresholds
display(['Clustering ',char(chan)])

finalVar1Thresh=mean(allSampleThresholds(:,1));
finalVar2Thresh=mean(allSampleThresholds(:,2));
OFF_clust_points=intersect(find(clusterVar1<finalVar1Thresh),find(clusterVar2<finalVar2Thresh));
ON_clust_points=setdiff(1:1:length(clusterVar1),OFF_clust_points)';
clusterIDX_total=ones(length(vsPNE),1);
clusterIDX_total(ON_clust_points)=2;

%%%%% Find all OFF periods
OFFgaps=find(diff(vsPNEtime(OFF_clust_points))>(1/floor(OFFDATA_var.PNEfs))); %find last epoch of each episode
numOFF=length(OFFgaps); %number of OFF periods
OFFperiod=[];
%loop going through all OFF periods
for ep = 1:numOFF
    if ep==1
        StartOFF=OFF_clust_points(1); %find start point of this OFF period
        EndOFF=OFF_clust_points(OFFgaps(ep)); %find last point of this OFF period
    elseif ep==numOFF
        StartOFF=OFF_clust_points(OFFgaps(ep)+1); %find start point of this OFF period
        EndOFF=OFF_clust_points(end); %find last point of this OFF period
    else
        StartOFF=OFF_clust_points(OFFgaps(ep-1)+1); %find start point of this OFF period
        EndOFF=OFF_clust_points(OFFgaps(ep)); %find last point of this OFF period
    end
    OFFperiod=[OFFperiod; vsPNEtime(StartOFF) vsPNEtime(EndOFF)];
    clear startOFF endOFF
end
OFFDATA.(chan).(stage)=OFFperiod;

if strcmp(fieldnames(OFFDATA.(chan)),'AllOFFtimes')==0
    OFFDATA.(chan).AllOFFtimes=[];
end
try
    OFFDATA.(chan).AllOFFtimes=[OFFDATA.(chan).AllOFFtimes,OFF_clust_points];
end

clear OFFperiod OFF_clust_points ON_clust_points clusterVar1 clusterVar2 ...
      clusterIDX_total clustVar1ThreshScaled clustVar2ThreshScaled

end
end

[OFFDATA]=OFFAD_channelstats(OFFDATA,OFFDATA_var);


%% Function used to optimise thresholds
function falsePos=clusterThresholds(thresholds)
global X
global idx
global OFFclusterINDEX
%%% Compare threshold cluster with original cluster
cluster_one_sample=intersect(find(X(:,1)<thresholds(1)),find(X(:,2)<thresholds(2)));
cluster_two_sample=setdiff(1:1:length(find(X(:,1))),cluster_one_sample)';
clusterIDX_sample=ones(length(X(:,1)),1);
clusterIDX_sample(cluster_two_sample)=2;
%set original clustering index to 2 clusters only (combine medium and large clusters)
IDX_2clust=ones(length(idx),1);
IDX_2clust(idx~=OFFclusterINDEX)=2;
falsePos=sum(abs(diff([clusterIDX_sample,IDX_2clust],1,2)));




