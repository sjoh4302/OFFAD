function [OFFDATA]=OFFAD_clustering(importDataVar)
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
OFFDATA.selectStages=[importDataVar(24).Value,...      
              importDataVar(26).Value,...
              importDataVar(28).Value];
OFFDATA.ignoreChannels=double(string(regexp(importDataVar(30).String,'\d*','match')));
OFFDATA.clustVar1Select=double(string(importDataVar(32).Value)); %1 = amplitude, 2=power
OFFDATA.clustVar1Smooth=double(string(importDataVar(34).String));
OFFDATA.clustVar2Select=double(string(importDataVar(35).Value)); %1 = amplitude, 2=power
OFFDATA.clustVar2Smooth=double(string(importDataVar(37).String));
OFFDATA.percSamp=double(string(importDataVar(39).String)); %percentage of signal to sample for cluster thresholding

%Close import window
g.Import = findobj('tag', 'OFFAD_IMPORT');
close(g.Import)


%%%%% Perform clustering
%Choose stage to analyse
allstages=["nr","r","w"];
stages=allstages(logical(OFFDATA.selectStages));
for stageNum=1:length(stages)
stage=stages(stageNum);
%% collect bin number for start and end of each NREM epoch
% cutting off one 4 second epoch at the beginning and end of each episode to exclude transitional states

%load vigilance state information for this animal and recording date ('nr' are all artefact free NREM epochs)
vigState = load(OFFDATA.VSpathin,'-mat',stage);
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

exampleObject = matfile(OFFDATA.PNEpathin);
channelNums = who(exampleObject);
%Sort channels
numOnly=double(string(regexp(string(channelNums),'\d*','match')));
[chanOrderNums,chanOrder]=sort(numOnly);
%Remove selected channels
if ~isempty(OFFDATA.ignoreChannels)
    chanOrder=(chanOrder(~abs(sum(chanOrderNums==OFFDATA.ignoreChannels,2))));
end

PNElength=size(exampleObject,char(channelNums(1)),2);
OFFDATA.(stage).StartOP=sparse(repmat(logical(0),PNElength,length(chanOrder)));
OFFDATA.(stage).EndOP=sparse(repmat(logical(0),PNElength,length(chanOrder)));
OFFDATA.(stage).AllOP=sparse(repmat(logical(0),PNElength,length(chanOrder)));
OFFDATA.ChannelsFullName=string(channelNums(chanOrder));
OFFDATA.Channels=double(string(cellfun(@(X) regexp(X,'\d*','match'),channelNums(chanOrder),'UniformOutput',false)));
     
clear exampleObject PNElength

for chanNum = 1:length(chanOrder)
    chan=string(channelNums(chanOrder(chanNum)));    
    
    %skip channel if it's not a good channel (was excluded on visual inspection due to bad signals)
            
        %%% load pNe signal
        PNE = load(OFFDATA.PNEpathin,'-mat',chan);
        PNE = PNE.(chan);
        
        PNE = abs(PNE); %take absolute values
        
        if OFFDATA.PNEunit=="V" %Convert to uV
            PNE=PNE/1000000;
        end
        PNE = PNE(1:1000000);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%REMOV!!!!!!!!!!!!!!!!!
                
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

%% Spectral clustering (3 clusters)
%Extracts 30s samples of recording and perfroms spectral clustering
sampleLength=30; %in seconds
totalSamp=floor(((length(vsPNE)/OFFDATA.PNEfs)/sampleLength));
sampSize=round(OFFDATA.PNEfs*sampleLength);
sampStarts=([0:1:totalSamp-1]*sampSize)+1;
numSamp=round(totalSamp*(OFFDATA.percSamp/100));
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
OFFgaps=find(diff(vsPNEtime(OFF_clust_points))>1); %find last epoch of each episode
numOFF=length(OFFgaps); %number of OFF periods
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
    StartOFF=vsPNEtime(StartOFF);
    EndOFF=vsPNEtime(EndOFF);
    OFFDATA.(stage).StartOP(StartOFF,chanNum)=1;
    OFFDATA.(stage).EndOP(EndOFF,chanNum)=1;
    OFFDATA.(stage).AllOP(StartOFF:EndOFF,chanNum)=1;
   
end

clear OFFperiod OFF_clust_points ON_clust_points clusterVar1 clusterVar2 ...
      clusterIDX_total clustVar1ThreshScaled clustVar2ThreshScaled

end
end

[OFFDATA]=OFFAD_channelstats(OFFDATA);


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



%Draw a heatmap 
figure
  
    Var1Max=round(prctile(clusterVar1,99));
    Var2Max=round(prctile(clusterVar2,99));
    nhistBins=800;
    Var1edges=round(min(clusterVar1)):Var1Max/nhistBins:Var1Max;
    Var2edges=round(min(clusterVar2)):Var2Max/nhistBins:Var2Max;
[N,Var1edges,Var2edges] = histcounts2(clusterVar1,clusterVar2,Var1edges,Var2edges);
   
    imagesc(flip(N'))
    oldX=xticks;
    xticklabels(oldX/max(oldX)*ampMax)
    oldY=yticks;
    yticklabels(flip(oldY/max(oldY)*scaloMax))
    c = hot;
    c = flipud(c);
    colormap(c);
    %colormap('hot','Direction','reverse')
    %caxis([50 430758])
    %set(gca,'ColorScale','log')
    %xlim([-0.05 1.05])
    
    randStart=randi(length(clusterVar1-100000),1);
    %sampCluster=[clusterVar1(randStart:randStart+100000)',clusterVar2(randStart:randStart+100000)'];
    randPoints=randi(length(clusterVar1),100000,1);
    sampCluster=[clusterVar1(randPoints)',clusterVar2(randPoints)'];
    sampCluster=[rescale(sampCluster(:,1)),rescale(sampCluster(:,2))];
    clusterData=[clusterVar1',clusterVar2'];
    clusterData=[rescale(clusterData(:,1)),rescale(clusterData(:,2))];
    
    randPoints=randi(length(clusterVar1_full),100000,1);
    sampCluster=[clusterVar1_full(randPoints)',clusterVar2_full(randPoints)'];
    clusterData=[clusterVar1_full',clusterVar2_full'];
    %%
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
    eva = evalclusters(sampCluster,allIDX,'CalinskiHarabasz')
    %eva = evalclusters(sampCluster,'gmdistribution','CalinskiHarabasz','KList',[1:5])
    % eva2 = evalclusters(clusterData2,allIDX2,'CalinskiHarabasz')
    randPoints=randi(length(clusterVar1),1000000,1);
    sampCluster=[clusterVar1(randPoints)',clusterVar2(randPoints)'];
    options = statset('MaxIter',1000,'TolFun',1e-5);
    GMModel = fitgmdist(sampCluster,eva.OptimalK,'Replicates',20,'Options',options);
    IDX = cluster(GMModel,clusterData);
    
    [~,offCluster]=min(GMModel.mu(:,1));
    startPoint=round(732*OFFDATA.PNEfs);
    figure
    stem(PNE(startPoint:startPoint+4000),'-','LineWidth',1,'Marker','none','ShowBaselin','off')
    hold on
    stem(find(IDX(startPoint:startPoint+4000)==offCluster),PNE(startPoint-1+find(IDX(startPoint:startPoint+4000)==offCluster)),'-','LineWidth',1,'Marker','none','ShowBaselin','off','Color','red','LineStyle','-')
      %%  
    %eva = evalclusters([clusterVar1(sampClusts)',clusterVar2(sampClusts)'],...
     %   'gmdistribution','CalinskiHarabasz','KList',[1:4])
    
    
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
    figure('Position',[80,150,1400,600])
    subplot(1,2,2)
    d = 700; % Grid density
    x1 = linspace(min(clusterVar1), max(clusterVar1), d);
    x2 = linspace(min(clusterVar2), max(clusterVar2), d);
    [x1grid,x2grid] = meshgrid(x1,x2);
    X0 = [x1grid(:) x2grid(:)];
    mahalDist = mahal(GMModel,X0);
    threshold = sqrt(chi2inv(0.99,2));
    hold on
    [~,clustOrder]=sort(GMModel.mu(:,1)); 
    for m = 1:size(mahalDist,2)
        currentClust=clustOrder(m);
        idx = mahalDist(:,currentClust)<=threshold;
        h2 = plot(X0(idx,1),X0(idx,2),'.','MarkerSize',1);
        uistack(h2,'bottom');
    end 
    
    %P = posterior(GMModel,a);
     
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
   % caxis([0 2000])
    %set(gca,'ColorScale','log')
    %set(gca,'CLim',[0.1 2000])
    
    figure
    dispSample=randi(length(clusterData),1,10000);
    gscatter(clusterData(dispSample,1),clusterData(dispSample,2),IDX(dispSample))
    xlim([clusterAxes.XLim(1),clusterAxes.XLim(2)])
    ylim([clusterAxes.YLim(1),clusterAxes.YLim(2)])