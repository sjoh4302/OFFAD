% 
clear all
close all
tic
%% Intialise variables
for x=1:4
if x==1
%%% Recording information 1
provider='Lukas';%'Mathilde';
derivation='NA';%'fro';
mousename='MH';%'Ph';
BLdate='20190217';%'200318';
LightPhase='DL';%'L' %(light period)
AllChannels=[1:8];
GoodChannels=[1:8]; %select which channels to analyse
fs=498.2462; %sampling rate of the pNe signal
epochlength=4; %vigilance state scoring in 4s epochs
elseif x==2
%%% Recording information 2
provider='Lukas';%'Mathilde';
derivation='NA';%'fro';
mousename='MH';%'Ph';
BLdate='20190217';%'200318';
LightPhase='DL';%'L' %(light period)
AllChannels=[9:16];
GoodChannels=[9:16]; %select which channels to analyse
fs=498.2462; %sampling rate of the pNe signal
epochlength=4; %vigilance state scoring in 4s epochs
elseif x==3
%%% Recording information 3
provider='Mathilde';
derivation='fro';
mousename='Ph';
BLdate='200318';
LightPhase='L' %(light period)
AllChannels=[1:6];
GoodChannels=[1:6]; %select which channels to analyse
fs=498.2462; %sampling rate of the pNe signal
epochlength=4; %vigilance state scoring in 4s epochs
elseif x==4
%%% Recording information 4
provider='Mathilde';
derivation='fro';
mousename='Qu';
BLdate='200318';
LightPhase='L' %(light period)
AllChannels=[3:8];
GoodChannels=[3:8]; %select which channels to analyse
fs=498.2462; %sampling rate of the pNe signal
epochlength=4; %vigilance state scoring in 4s epochs
end

%%% File and path description
path='D:\DPhil\Off-period detection\';
pathinPNE=[path,'Off-period detection\']; %path for pNe files
pathinVS=[path,'TestSignals\']; %path for VS files

%% collect bin number for start and end of each NREM epoch
% cutting off one 4 second epoch at the beginning and end of each episode to exclude transitional states

%load vigilance state information for this animal and recording date ('nr' are all artefact free NREM epochs)
filenameVS=['Data_',provider,'\',mousename,'_',BLdate,'_',LightPhase,'_',derivation,'_VSspec'];
load([pathinVS,filenameVS,'.mat'],'-mat','nr');

% find NREM episodes
endEpi=find(diff(nr)>1); %find last epoch of each episode
endEpi=[endEpi;length(nr)]; %add final episode end
startEpi=find(diff(nr)>1)+1; %find first epoch of each episode
startEpi=[1;startEpi]; %add first episode start
numepi=length(startEpi); %number of NREM episodes
cleanepochs=[];

%loop going through all NREM episodes 
for ep = 1:numepi
    
    Startepoch=nr(startEpi(ep)); %find start epoch of this NREM episode
    Endepoch=nr(endEpi(ep)); %find second last epoch of this NREM episode
    
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

for chan = AllChannels(1):AllChannels(end)
    
    %skip channel if it's not a good channel (was excluded on visual inspection due to bad signals)
    if ~ismember(chan,GoodChannels)
        continue
    else
        
        %%% load pNe signal
        filename=['Data_',provider,'\','pNe_data_',mousename,'_',BLdate,'_',LightPhase,'_Channels',num2str(AllChannels(1)),'to',num2str(AllChannels(end))];
        sig=load([pathinVS,filename,'.mat'],'-mat',['Ch',num2str(chan)]);
        sig = sig.(['Ch',num2str(chan)]);
        
        %sig=sig/1000000; %transform into uV values
        
        sig=abs(sig); %take absolute values
          
        
        
        
        %%% collect pNe signal for all NREM episodes
        NremSig=NaN(1,length(sig)); %NaN vectors to be filled with NREM pNe signal
        
        %loop going through all NREM epochs (except for first and last)
        for ep = 1:numepochs
            
            Startepoch=cleanepochs(ep,1); %start of respective epoch
            Endepoch=cleanepochs(ep,2); %end of respective epoch
            
            %fill NaN vectors with pNe signal for this NREM episode
            StartBin=ceil(Startepoch*epochlength*fs);
            EndBin=floor((Endepoch-1)*epochlength*fs);
            NremSig(StartBin:EndBin)=sig(StartBin:EndBin);
            
        end
        %concatenate the signal from all selected NREM episodes to get rid of gaps
        NremSig=NremSig(~isnan(NremSig));
        
        
    end
%end



%% Data processing
%%% Smooth pNe signal
% Convolve absolute NREM signal with a gaussian window 
winSize=30; %in data points 5 for Lukas
smoothNremSig=conv(NremSig,gausswin(winSize))';
smoothNremSig=smoothNremSig(15:end-15);

%%% Generate scalogram of raw pNe signal
percOverlap=10; %overlap between segments
LB_freq=40; %find signal power from 0-LB-freq
[scalogramWT,F] = offPeriodScalogram(NremSig,fs,percOverlap,LB_freq);


%% Spectral clustering (3 clusters)
%tic
% Select plotting option
plotting=0; %0 = dont plot, 1 = plot

%Extracts a 30s sample of NREM recording and perfroms spectral clustering
sampleLength=30; %in seconds
percSamp=10; %percentage of signal to sample for cluster thresholding
numSamp=round(((length(NremSig)/fs)/sampleLength)*(percSamp/100)); 
color = lines(3);
sampSize=round(fs*sampleLength);

%Declare global variables for thresholding function
global X
global idx
global OFFclusterINDEX

for n = 1:numSamp 
startSeed=randi(length(smoothNremSig)-sampSize,1);

X=[smoothNremSig(startSeed:startSeed+sampSize-1),scalogramWT(startSeed:startSeed+sampSize-1)'];
Xscaled=[rescale(X(:,1)),rescale(X(:,2))];

idx = spectralcluster(Xscaled,3);

if plotting==1
    
    figure(n)
    subplot(3,1,1)
    gridx1=0:0.01:1;
    gridx2=0:0.01:1;
    [x1,x2] = meshgrid(gridx1, gridx2);
    x1 = x1(:);
    x2 = x2(:);
    xi = [x1 x2];
    %ksdensity(Xscaled,'PlotFcn','surf');
    ampMax=1;
    scaloMax=1;
    nhistBins=100;
    N=hist3(Xscaled,...
        'Ctrs',{0:ampMax/nhistBins:ampMax 0:scaloMax/nhistBins:scaloMax},...
        'CDataMode','auto','FaceColor','interp');
    %N=hist3([smoothNremSig,scalogramWT(1:length(NremSig))'],[50 50]);
    %subplot(3,1,3)
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
    set(gca,'ColorScale','log')
    %xlim([-0.05 1.05])


    subplot(3,1,2)
    %for i = 1:16
        %idx = dbscan(Xscaled,0.0145,100);
        %idx = dbscan(Xscaled,0.036,1000);
        %idx = spectralcluster(Xscaled,3); %good at startSeed=12152610
        %idx = kmeans(Xscaled(:,1),3);
       % subplot(4,4,i)
        gscatter(Xscaled(:,1),Xscaled(:,2),idx,color);
    %end

    subplot(3,1,3)
    displayLength=5000;
    currentSegment=NremSig(startSeed:startSeed+displayLength-1);
    sampleIDX=idx(1:displayLength);
    for  j = 1:max(unique(sampleIDX))
    %opCluster=2j;
        OPtimes=find(sampleIDX==j);
        OPgaps=find(diff(OPtimes)>1); %find last epoch of each episode
        numOP=length(OPgaps); %number of Off-period episodes (-1 because last epoch doesn't end with a gap)
        off_periods=[];
        for ep = 1:numOP
            if ep==1
                StartTime=OPtimes(ep);
                EndTime=OPtimes(OPgaps(ep));
            else
                StartTime=OPtimes(OPgaps(ep-1)+1); %find start epoch of this Off-period 
                EndTime=OPtimes(OPgaps(ep)); %find second last epoch of this Off-periode
            end
            off_periods=[off_periods; StartTime EndTime];

        end

        for i=1:length(off_periods)
            plot(off_periods(i,1):off_periods(i,2),currentSegment(off_periods(i,1):off_periods(i,2)),'color',color(j,:))
            hold on
        end
    end
    oldXlabels=xticks;
    xticklabels(round(oldXlabels/fs))
    
end

%%% Find thresholds from clusterings
[smoothThreshScaled,OFFclusterINDEX]=min([max(Xscaled(idx==1,1)),max(Xscaled(idx==2,1)),max(Xscaled(idx==3,1))]);
smoothThresh=X(find(Xscaled(:,1)==smoothThreshScaled),1);
[scalogramThreshScaled]=min([max(Xscaled(idx==1,2)),max(Xscaled(idx==2,2)),max(Xscaled(idx==3,2))]);
scalogramThresh=X(find(Xscaled(:,2)==scalogramThreshScaled),2);

% Use fminsearch to find thresholds which minimse the loss function (false
% positive plus false negative
x0 = [smoothThresh,scalogramThresh];
[sampleThresh,fval] = fminsearch(@clusterThresholds,x0);
allSampleThresholds(n,:)=sampleThresh;
end
%toc
%% Find off periods using clustered thresholds

finalSmoothThresh=mean(allSampleThresholds(:,1));
finalScalogramThresh=mean(allSampleThresholds(:,2));
cluster_one_total=intersect(find(smoothNremSig<finalSmoothThresh),find(scalogramWT<finalScalogramThresh));
cluster_two_total=setdiff(1:1:length(smoothNremSig),cluster_one_total)';
clusterIDX_total=ones(length(NremSig),1);
clusterIDX_total(cluster_two_total)=2;

%%%% Show example of final thresholding 
%figure(5)
%plot(cluster_one_total,NremSig(cluster_one_total),'.')
%hold on
%plot(cluster_two_total,NremSig(cluster_two_total),'.')

%Find all OFF periods
OFFgaps=find(diff(cluster_one_total)>1); %find last epoch of each episode
numOFF=length(OFFgaps); %number of OFF periods
OFFperiod=[];

%loop going through all OFF periods
for ep = 1:numOFF
    if ep==1
        StartOFF=cluster_one_total(1); %find start point of this OFF period
        EndOFF=cluster_one_total(OFFgaps(ep)); %find last point of this OFF period
    elseif ep==numOFF
        StartOFF=cluster_one_total(OFFgaps(ep)+1); %find start point of this OFF period
        EndOFF=cluster_one_total(end); %find last point of this OFF period
    else
        StartOFF=cluster_one_total(OFFgaps(ep-1)+1); %find start point of this OFF period
        EndOFF=cluster_one_total(OFFgaps(ep)); %find last point of this OFF period
    end
    OFFperiod=[OFFperiod; StartOFF EndOFF];
    clear startOFF endOFF
    
end
savename=[pathinVS,'Data_',provider,'\','OFFperiod_',mousename,'_',BLdate,'_',LightPhase,'_Channel',num2str(chan),'_1%'];       
save(savename,'OFFperiod','allSampleThresholds')

clearvars -except provider derivation mousename BLdate LightPhase ...
    AllChannels GoodChannels fs epochlength path pathinPNE pathinVS ...
    filenameVS gaps numepi cleanepochs Startepoch Endepoch nr numepochs
end
end

%% Analyse off periods
% fs2=round(fs);
% OFFperiod_durations=[OFFperiod(:,2)-OFFperiod(:,1)+1]/fs2*1000; %in milliseconds
% ONperiod_durations=[OFFperiod(2:end,1)-OFFperiod(1:end-1,2)]/fs2*1000; %in milliseconds
% numepochs=length(cleanepochs);
% meanDuration_unprocessed=mean(OFFperiod_durations)
% figure(10)
% binEdges=[0:(1/fs2)*1000*2:1000];
% subplot(1,2,1)
% histogram(OFFperiod_durations,binEdges)
% subplot(1,2,2)
% h=histogram(ONperiod_durations,binEdges); 
% splitMs=h.BinEdges(find(islocalmin(h.BinCounts,'MaxNumExtrema',1)))+h.BinWidth/2; %Point of bimodality
% ONbimodalSplit=round(splitMs/1000/(1/fs2));
% 
% figure(11)
% subplot(2,1,1)
% randPos=randi(length(OFFperiod)-20);
% numOFFplt=[randPos,randPos+20];
% plot(OFFperiod(numOFFplt(1),1):OFFperiod(numOFFplt(2),2),NremSig(OFFperiod(numOFFplt(1),1):OFFperiod(numOFFplt(2),2)))
% hold on
% for i=numOFFplt(1):numOFFplt(2)
%             plot(OFFperiod(i,1):OFFperiod(i,2),NremSig(OFFperiod(i,1):OFFperiod(i,2)),'color','r')
%             hold on
% end
% 
% %Plot OFF-periods with short interuptions removed
% subplot(2,1,2) 
% plot(OFFperiod(numOFFplt(1),1):OFFperiod(numOFFplt(2),2),NremSig(OFFperiod(numOFFplt(1),1):OFFperiod(numOFFplt(2),2)))
% hold on
% for i=numOFFplt(1):numOFFplt(2)
%             plot(OFFperiod(i,1):OFFperiod(i,2),NremSig(OFFperiod(i,1):OFFperiod(i,2)),'color','r')
%             hold on
%             if i<numOFFplt(2)
%                 if length(NremSig(OFFperiod(i,2):OFFperiod(i+1,1)))<ONbimodalSplit
%                    plot(OFFperiod(i,2):OFFperiod(i+1,1),NremSig(OFFperiod(i,2):OFFperiod(i+1,1)),'color','r')
%                 else
%                    plot(OFFperiod(i,2):OFFperiod(i+1,1),NremSig(OFFperiod(i,2):OFFperiod(i+1,1)),'color','b')
%                 end
%             end
% end
% 
% concatOFFperiod=OFFperiod(1,:);
% concatCounter=1;
% for i=2:length(OFFperiod)
%     if (OFFperiod(i,1)-concatOFFperiod(concatCounter,2))<ONbimodalSplit
%         concatOFFperiod(concatCounter,2)=OFFperiod(i,1);
%     else
%         concatOFFperiod=[concatOFFperiod;OFFperiod(i,:)];
%         concatCounter=concatCounter+1;
%     end
% end
% concatOFFperiod_duration=[concatOFFperiod(:,2)-concatOFFperiod(:,1)+1]/fs2*1000; %in milliseconds
% meanDuration_concat=mean(concatOFFperiod_duration)
% figure(12)
% binEdges=[0:(1/fs2)*1000*2:1000];
% subplot(1,2,1)
% histogram(OFFperiod_durations,binEdges)
% subplot(1,2,2)
% histogram(concatOFFperiod_duration,binEdges)

%%%Use code below to plot discontinous episodes, slow and plot hard to
%%%manipulate due to size of data
% for  j = 1
% 
%     ClustTimes=find(clusterIDX==j);
%     ClustGaps=find(diff(ClustTimes)>1); %find last epoch of each episode
%     numEp=length(ClustGaps); %number of cluster episodes (-1 because last epoch doesn't end with a gap)
%     currentCluster=[];
%     for ep = 1:numEp
%         if ep==1
%             StartTime=ClustTimes(ep);
%             EndTime=ClustTimes(ClustGaps(ep));
%         else
%             StartTime=ClustTimes(ClustGaps(ep-1)+1); %find start epoch of this cluster episode 
%             EndTime=ClustTimes(ClustGaps(ep)); %find second last epoch of this cluster episode
%         end
%         currentCluster=[currentCluster; StartTime EndTime];
% 
%     end
%     
%     for i=1:length(currentCluster)
%         plot(currentCluster(i,1):currentCluster(i,2),NremSig(currentCluster(i,1):currentCluster(i,2)),'color',color(j,:))
%         hold on
%     end
% end
% 




%plot(find(X(:,1)>smoothThresh & X(:,2)>scalogramThresh),NremSig(find(X(:,1)>smoothThresh & X(:,2)>scalogramThresh)+startSeed),'.')
%hold on
%plot(find(X(:,1)<smoothThresh & X(:,2)<scalogramThresh)+startSeed,NremSig(find(X(:,1)<smoothThresh & X(:,2)<scalogramThresh)+startSeed),'.')


%plot(find(X(:,1)>smoothThresh & X(:,2)>scalogramThresh)+startSeed,NremSig(find(X(:,1)>smoothThresh & X(:,2)>scalogramThresh)+startSeed))
%hold on
%plot(find(X(:,1)<smoothThresh & X(:,2)<scalogramThresh)+startSeed,NremSig(find(X(:,1)<smoothThresh & X(:,2)<scalogramThresh)+startSeed))
toc
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
end