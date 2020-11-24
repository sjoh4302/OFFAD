%%%%
clear all
close all
%{
GENERAL INFORMATION
This script aims to define an 'noise threshold' in each individual pNe channel based on bimodality the spike amplitude distribution for subsequent analysis of
OFF-periods

AUTHOR INFORMATION
- script written by LBK 20200706
    (contact:lukas.krone@dpag.ox.ac.uk)
- last update: ----

INPUT
- mousename
- baseline date 'BL date'
- 'GoodChannels' for each animal (see animalinfo Excel-sheet)
- pNe signal for each channel
- vigilance state info from EEG scoring
- start and end time (in s) of time window (ideally 10 s) for 'test signal'
with clear ON and OFF states during NREM sleep

OUTPUT
- figure showing thresholds from differnt approaches on test signal
%}

%% clear all variables, close all figures
%clear all
%close all

%% file and path description
path='D:\DPhil\Off-period detection\';
pathinPNE=[path,'Off-period detection\']; %path for pNe files
pathinVS=[path,'TestSignals\']; %path for VS files


%% recording information
mousename='Ph';
%mousename='MH';
BLdate='200318';
%BLdate='20190217';
LightPhase='L'; %(light period)
GoodChannels=[8:14];
%GoodChannels=[1:6 8:14];


SR=498.2462; %sampling rate of pNe signal
f1=SR;
epochlength=4; %vigilance state scoring in 4s epochs


%% provide time window of 10s NREM sleep as 'test signal' for visualisation (unit seconds, find good trace in OpenScope)
TestWindowStart=9675;  %start time in sec; %TestWindowStart=round(TestWindowStart*SR); find respective pNe bin
TestWindowEnd=9685;  %end time in sec; %TestWindowEnd=round(TestWindowEnd*SR); find respective pNe bin
TestWindowLength=TestWindowEnd-TestWindowStart; %length in sec
% TestWindow=[TestWindowStart:TestWindowEnd];

%% collect bin number for start and end of each NREM epoch
% cutting off one 4 second epoch at the beginning and end of each episode to exclude transitional states

%load vigilance state information for this animal and recording date ('nr' are all artefact free NREM epochs)
%filenameVS=[mousename,'-',BLdate,'-eeg17','-VSspec'];
filenameVS=['Data_Mathilde\',mousename,'_',BLdate,'_',LightPhase,'_fro_VSspec'];
load([pathinVS,filenameVS,'.mat'],'-mat','nr');

% find NREM episodes
gaps=find(diff(nr)>1); %find last epoch of each episode
numepi=length(gaps); %number of NREM episodes (-1 because last epoch doesn't end with a gap)
cleanepochs=[];

%loop going through all NREM episodes (exept for first and last)
for ep = 2:numepi
    
    Startepoch=nr(gaps(ep-1)+1); %find start epoch of this NREM episode
    Endepoch=nr(gaps(ep))-1; %find second last epoch of this NREM episode
    
    if Endepoch-Startepoch<3 %if NREM episode 2 epochs or shorter, there'll be no signal because w're cutting the first and last epoch
        continue
    else %if NREM episode at least 3 epochs, find bins corresponding to the start of the second epoch and end of the second last epoch
        cleanepochs=[cleanepochs; Startepoch Endepoch];
    end
end

numepochs=length(cleanepochs);

%% load and process pNe signal for each channel

%{
1) convert it into absolute values
2) create resampled copy with new sampling rate of 256 Hz
3) compare absolute signal before and after resampling
%}
SmoothedSpikeDistribution=NaN(16,151);


for chan = 1:16
    
    %skip channel if it's not a good channel (was excluded on visual inspection due to bad signals)
    if ~ismember(chan,GoodChannels)
        continue
    else
        
        %load pNe signal
        %filename=[mousename,'-pNe-',BLdate,'-ch',num2str(chan)];
        %load([pathinVS filename '.mat']);
        filename=['Data_Mathilde\','pNe_data_',mousename,'_',BLdate,'_',LightPhase,'_Channels1to6'];
        sig=load([pathinVS,filename,'.mat'],'-mat',['Ch',num2str(chan)])
        sig = sig.(['Ch',num2str(chan)]);
        
        %%%%%transform
        %sig=sig/1000000; %transform into uV values
        
        sig=abs(sig); %take absolute values
        
        
        %resample to 256 Hz
        % f1=SR;
        % f2=256;
        % num_chunks=4;
        % tail_length=20000;
        % visualize=0;
        %
        % sig_resamp=resampling(sig,f1,f2,num_chunks,tail_length,visualize);
        
        
        %plot figure to see resampling results
        indSRorig=round(TestWindowStart*f1):1:round(TestWindowEnd*f1);
        xOrig=10/length(indSRorig):10/length(indSRorig):TestWindowLength;
        testsignal=sig(indSRorig);
        
        % indSRresamp=round(TestWindowStart*f2):1:round(TestWindowEnd*f2);
        % xResamp=10/length(indSRresamp):10/length(indSRresamp):TestWindowLength;
        % testsignal_resamp=sig_resamp(indSRresamp);
        
        % plot testsignal as figure(1)
        figure(1)
       % subplot(4,4,chan)
        p1=plot(xOrig,testsignal,'b');
        % lgd_p1='original signal';
        hold on
        % p2=plot(xResamp,testsignal_resamp,'r');
        % lgd_p2='resampled testsignal';
        
        %figure settings
        % legend(lgd_p1,lgd_p2)
        ylabel('Signal amplitude (in uV)')
        xlabel('Time [s]')
        title(['test signal - channel: ',num2str(chan)])
        ylim([0 200])
        % set(gcf,'Position',[200 200 1500 500])
        
        
        
        %% collect pNe signal for all NREM episodes
        %NaN vectors to be filled with NREM pNe signal
        NremSig=NaN(1,length(sig));
        
        %loop going through all NREM epochs (exept for first and last)
        for ep = 1:numepochs
            
            Startepoch=cleanepochs(ep,1); %start of respective epoch
            Endepoch=cleanepochs(ep,2); %end of respective epoch
            
            %fill NaN vectors with pNe signal for this NREM episode
            StartBin=ceil(Startepoch*epochlength*f1);
            EndBin=floor(Endepoch*epochlength*f1);
            NremSig(StartBin:EndBin)=sig(StartBin:EndBin);
            
        end
        %concatenate the signal from all selected NREM episodes to get rid of gaps
        totsig=NremSig(~isnan(NremSig));
        
        
        %% threshold determination for respective channel - try to identify optimal cutoff in bimodal distribution of spike amplitudes or in smoothed test signal using different smoothing factors
        %{
    1) plot histogram of spike amplitudes
    2) find cutoff point in bimodal distribution
    3) set threshold manually
        %}
        
        for av=1:16
            
            %averagetotsig=movmean(totsig,av); %smooth total signal
            averagetotsig=movmean(totsig,av); %smooth total signal
            %averagetotsig=conv(totsig,gausswin(av)); %smooth total signal with Gaussian
            
            
            figure(101+chan)
            subplot(4,4,av)
            h=histc(averagetotsig,[0:1:150]);
            bar(h);
            %     set(gca, 'XTick', [0:5:30]);
            %     xticklabels({[0:25:150]});
            xlabel('Spike Amplitude [uV]')
            
            if av==5
                SmoothedSpikeDistribution(chan,:)=h;
            end
            
            %for visualisation only
            %{
    averagetestsig=movmean(testsignal,av); %smooth test signal
    figure(11)
    subplot(4,4,av)
    plot(averagetestsig,'b')
            %}
        end
        sgtitle(['Histogram of spike amplitudes - channel: ',num2str(chan)])
        %suptitle({['distribution of spike amplitutes in NREM signal'],['with increasing moving window size']})
        
        
        %%%% activate in first round to find minimum in bimodal distribution %%%%
        % (later import vector from folder)
                prompt1 = 'Threshold based on inspection of spike distribution (bimodality):';
                threshold_1(chan) = input(prompt1);
        
        %plot visual threshold
        figure(1)
        %hold on
        p=plot([0 TestWindowLength],[threshold_1(chan) threshold_1(chan)],'g--','LineWidth',2);
        %         lgd_p4='threshold based on inspection of spike distribution';
        %         legend(lgd_p1,lgd_p2,lgd_p3,lgd_p4)
        
        
        
    end
end

%% Clustering

%Find spectral power of concatenated NREM signal using scalograms



testSig6=sig(indSRorig)';
testSig6Conv=conv(sig(indSRorig),gausswin(30))';

totsigConv=conv(totsig,gausswin(30))';
totsigConv=totsigConv(15:end-15);
clusters=kmeans(testSig6,3);
figure(89)
subplot(2,1,1)
scatter(testSig6(find(clusters==1)),ones(1,sum(clusters==1)))
hold on
scatter(testSig6(find(clusters==2)),ones(1,sum(clusters==2)))
scatter(testSig6(find(clusters==3)),ones(1,sum(clusters==3)))
subplot(2,1,2)
clusterThresh=min([max(testSig6(find(clusters==1))),max(testSig6(find(clusters==2))),...
    max(testSig6(find(clusters==3)))])
plot(testSig6)
hold on
p=plot([0 5000],[clusterThresh clusterThresh],'g--','LineWidth',2);

figure(31)
scatter(testSig6Conv(15:end-15),scalogramPower,'.')

numClust=3;
scalogramPower=mean(abs(WT(F>40,:)))';
clusters=kmedoids([testSig6Conv(15:end-15),scalogramPower],numClust,'Distance','sqeuclidean');
figure(90)
subplot(2,1,1)
for i =1:numClust
scatter(testSig6Conv((find(clusters==i))+14),scalogramPower(find(clusters==i)),'.')
hold on
clustPowerMean(i)=mean(scalogramPower(find(clusters==i)));
end
lowClust=find(clustPowerMean==min(clustPowerMean));
subplot(2,1,2)
scatter(1:length(testSig6Conv),testSig6Conv,'.')
hold on
scatter(find(clusters==lowClust)+14,testSig6Conv(find(clusters==lowClust)+14),'.')



X=[rescale(testSig6Conv(15:end-15)),rescale(scalogramPower)];
Z = linkage(X,'ward','euclidean');
figure(76)
subplot(3,1,1)
dendrogram(Z)
T = cluster(Z,'MaxClust',3);
subplot(3,1,2)
gscatter(X(:,1),X(:,2),T)
subplot(3,1,3)
scatter(1:length(X(:,1)),X(:,1),'.')
hold on
scatter(find(T==1),X(T==1,1),'.')

%Find off-period sizes
OPtimes=find(T==1);
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
figure(101)
subplot(2,1,1)
histogram(diff(off_periods,1,2),25)
subplot(2,1,2)
plot(testSig6)
hold on 
for i=1:length(off_periods)
    if diff(off_periods(i,:))>50
        plot(off_periods(i,1):off_periods(i,2),testSig6(off_periods(i,1):off_periods(i,2)),'r')
    end
end




options = statset('Display','final'); 
X=[rescale(testSig6Conv(15:end-15)),rescale(scalogramPower)];
gm = fitgmdist(X,2,'Options',options)

figure(98)
idx = cluster(gm,X);
cluster1 = (idx == 1); % |1| for cluster 1 membership
cluster2 = (idx == 2); % |2| for cluster 2 membership

gscatter(X(:,1),X(:,2),idx,'rb','+o')
legend('Cluster 1','Cluster 2','Location','best')



sig6Conv=conv(totsig,gausswin(30))';

scalogramPower=mean(abs(WT(F>16,:)))';
figure(91)
subplot(2,2,1)
histogram(sig6Conv,floor(max(sig6Conv)))
subplot(2,2,2)
plot(testSig6Conv)
hold on
plot([0 5000],[291 291],'g--','LineWidth',2);
subplot(2,2,3)
histogram(scalogramPower,150)
subplot(2,2,4)
plot(scalogramPower)
hold on
plot([0 5000],[2.6 2.6],'g--','LineWidth',2);


figure(15)
subplot(2,1,1)
scatter(totsigConv(:),scalogramWT(1:length(totsig)),'.')
subplot(2,1,2)
ampMax=round(max(totsigConv),-2);
scaloMax=round(max(scalogramWT),-1);
nhistBins=100;
N=hist3([totsigConv,scalogramWT(1:length(totsig))'],...
    'Ctrs',{0:ampMax/nhistBins:ampMax 0:scaloMax/nhistBins:scaloMax},...
    'CDataMode','auto','FaceColor','interp');
%N=hist3([totsigConv,scalogramWT(1:length(totsig))'],[50 50]);
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
caxis([50 430758])
set(gca,'ColorScale','log')


figure(104)
sampSize=10000;
%startSeed=randi(14572713-sampSize,1);
X=[totsigConv(startSeed:startSeed+sampSize-1),scalogramWT(startSeed:startSeed+sampSize-1)'];
subplot(3,1,1)
ampMax=round(max(X(:,1)),-2);
scaloMax=round(max(X(:,2)),-1);
nhistBins=100;
N=hist3(X,...
    'Ctrs',{0:ampMax/nhistBins:ampMax 0:scaloMax/nhistBins:scaloMax},...
    'CDataMode','auto','FaceColor','interp');
imagesc(flip(N'))
oldX=xticks;
xticklabels(oldX/max(oldX)*ampMax)
oldY=yticks;
yticklabels(flip(oldY/max(oldY)*scaloMax))
c = hot;
c = flipud(c);
colormap(c);

subplot(3,1,2)
idxCounter=1;
for e=1:100
    for mini=1:100
        idx= dbscan(X,e,mini);
        if max(idx)==2
            goodIDX(:,idxCounter)=idx;
            goodIDXsettings(idxCounter,:)=[e mini];
            idxCounter=idxCounter+1;
        end
    end
end
idx = dbscan(X,2,36);
%idx = spectralcluster(X,2);
numCluster=max(idx)+1;
color = lines(idx+2);
gscatter(X(:,1),X(:,2),idx,color);
xlim([0,ampMax])
%%
%%% Good example at startSeed=3513408, epsilon=0.04,min=500. bimodality in
%%% k-nearest n (500) at 0.02
%sampSize=10000; e=0.56,min=1000,2246099/3513408/10829563/553876/1815123 (~20% success rate)
%sampSize=20000;e=0.56,min=2000,(~50% success rate)
tic
figure(3)
sampSize=25000;
startSeed=randi(length(totsigConv)-sampSize,1);
X=[totsigConv(startSeed:startSeed+sampSize-1),scalogramWT(startSeed:startSeed+sampSize-1)'];
Xscaled=[rescale(X(:,1)),rescale(X(:,2))];
%subplot(4,1,1)
[nbors nborsD] = knnsearch(Xscaled,Xscaled,'k',1001);
meanDist=mean(nborsD(:,2:end),2);
maxDist=max(nborsD(:,2:end),[],2);
%loglog(sort(maxDist,'ascend'))
%xlim([sampSize/10 sampSize])
%ylim([0.001,1])
[counts,centr]=hist(meanDist,500);
%hist(meanDist,500)
[TF,P] = islocalmin(counts);
minima=centr(P==max(P));

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
%N=hist3([totsigConv,scalogramWT(1:length(totsig))'],[50 50]);
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
    idx = spectralcluster(Xscaled,3); %good at startSeed=12152610
    %idx = kmeans(Xscaled(:,1),3);
   % subplot(4,4,i)
    gscatter(Xscaled(:,1),Xscaled(:,2),idx,color);
%end

subplot(3,1,3)
displayLength=5000;
currentSegment=totsig(startSeed:startSeed+displayLength-1);
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
xticklabels(round(oldXlabels/f1))
toc

%%
noise_threshold=nanmean(threshold_1);
AveSmoothedSpikeDistributioin=nanmean(SmoothedSpikeDistribution);
[~,AveNoiselevel]=max(AveSmoothedSpikeDistributioin);

figure(99)
plot(SmoothedSpikeDistribution')
hold on
plot(nanmean(SmoothedSpikeDistribution),'k:','LineWidth',3)


for chan = 1:16
    
    %skip channel if it's not a good channel (was excluded on visual inspection due to bad signals)
    if ~ismember(chan,GoodChannels)
        continue
    else
        
        figure(1)
        subplot(4,4,chan)
        hold on
        t=plot([0 TestWindowLength],[noise_threshold noise_threshold],'r:','LineWidth',2);
        hold on
        n=plot([0 TestWindowLength],[AveNoiselevel+10 AveNoiselevel+10],'y:','LineWidth',2);
        
    end
end