close all
clear all

filename=('D:\DPhil\Off-period detection\TestSignals\Data_Lukas\OFFperiod_MH_20190217_DL_Channel');
for i = 1:16
    OFF_periods_Lukas(i)=load([filename,num2str(i),'_1%']);
end

filename=('D:\DPhil\Off-period detection\TestSignals\Data_Mathilde\OFFperiod_Ph_200318_L_Channel');
for i =1:6
     OFF_periods_MathildePh(i)=load([filename,num2str(i),'_1%']);
end

filename=('D:\DPhil\Off-period detection\TestSignals\Data_Mathilde\OFFperiod_Qu_200318_L_Channel');
for i =3:8
     OFF_periods_MathildeQu(i-2)=load([filename,num2str(i),'_1%']);
end

%% Plot ON-OFF histograms
%Matilde Ph
fs2=498;
figure
for i = 1:6
    OFFperiod=OFF_periods_MathildePh(i).OFFperiod;
    OFFperiod_durations=[OFFperiod(:,2)-OFFperiod(:,1)+1]/fs2*1000; %in milliseconds
    ONperiod_durations=[OFFperiod(2:end,1)-OFFperiod(1:end-1,2)]/fs2*1000; %in milliseconds
    binEdges=[0:(1/fs2)*1000*2:1000];
    subplot(2,6,i*2-1)
    histogram(OFFperiod_durations,binEdges)
    title(['Ph-Ch',num2str(i),' - OFF'])
    subplot(2,6,i*2)
    h=histogram(ONperiod_durations,binEdges); 
    title(['Ph-Ch',num2str(i),' - ON'])
end 

%Matilde Qu
fs2=498;
figure
for i = 1:6
    OFFperiod=OFF_periods_MathildeQu(i).OFFperiod;
    OFFperiod_durations=[OFFperiod(:,2)-OFFperiod(:,1)+1]/fs2*1000; %in milliseconds
    ONperiod_durations=[OFFperiod(2:end,1)-OFFperiod(1:end-1,2)]/fs2*1000; %in milliseconds
    binEdges=[0:(1/fs2)*1000*2:1000];
    subplot(2,6,i*2-1)
    histogram(OFFperiod_durations,binEdges)
    title(['Qu-Ch',num2str(i+2),' - OFF'])
    subplot(2,6,i*2)
    h=histogram(ONperiod_durations,binEdges); 
    title(['Qu-Ch',num2str(i+2),' - ON'])
end

%Lukas MH 1:8
fs2=498;
figure
for i = 1:8
    OFFperiod=OFF_periods_Lukas(i).OFFperiod;
    OFFperiod_durations=[OFFperiod(:,2)-OFFperiod(:,1)+1]/fs2*1000; %in milliseconds
    ONperiod_durations=[OFFperiod(2:end,1)-OFFperiod(1:end-1,2)]/fs2*1000; %in milliseconds
    binEdges=[0:(1/fs2)*1000*2:1000];
    subplot(2,8,i*2-1)
    histogram(OFFperiod_durations,binEdges)
    title(['MH-Ch',num2str(i),' - OFF'])
    subplot(2,8,i*2)
    h=histogram(ONperiod_durations,binEdges); 
    title(['MH-Ch',num2str(i),' - ON'])
end 


%Lukas MH 9:16
fs2=498;
figure
for i = 1:8
    OFFperiod=OFF_periods_Lukas(i+8).OFFperiod;
    OFFperiod_durations=[OFFperiod(:,2)-OFFperiod(:,1)+1]/fs2*1000; %in milliseconds
    ONperiod_durations=[OFFperiod(2:end,1)-OFFperiod(1:end-1,2)]/fs2*1000; %in milliseconds
    binEdges=[0:(1/fs2)*1000*2:1000];
    subplot(2,8,i*2-1)
    histogram(OFFperiod_durations,binEdges)
    title(['MH-Ch',num2str(i+8),' - OFF'])
    subplot(2,8,i*2)
    h=histogram(ONperiod_durations,binEdges); 
    title(['MH-Ch',num2str(i+8),' - ON'])
end

%Lukas MH 1:16
fs2=498;
figure
for i = 1:16
    OFFperiod=OFF_periods_Lukas(i).OFFperiod;
     ONperiod_durations=[OFFperiod(2:end,1)-OFFperiod(1:end-1,2)]/fs2*1000; %in milliseconds
    binEdges=[0:(1/fs2)*1000*2:1000];
    subplot(2,8,i)
    histogram(ONperiod_durations,binEdges); 
    title(['MH-Ch',num2str(i),' - ON'])
    
end

% %%
% for i = 1:16
%     MHthreshStats(i,1)=std(OFF_periods_Lukas(i).allSampleThresholds(:,1));
%     MHthreshStats(i,2)=std(OFF_periods_Lukas(i).allSampleThresholds(:,2));
%     MHthreshStats(i,3)=round(sum(diff(OFF_periods_Lukas(i).OFFperiod,1,2))/fs2);
% end
% % Bimodal ON channels in MH have more OFF periods 
% 
% 
% figure
% for i = 1:6
%     OFFperiod=OFF_periods_Lukas(i).OFFperiod;
%     for j=1:length(OFFperiod)
%         scatter(ones(diff(OFFperiod(i,:))),OFFperiod(i,1):OFFperiod(i,2))
%     end
% end

%% LFP analysis
fs_lfp=256;
%fs_lfp=305.1758;
fs_pne=498.246185;
epochlength=4;

%Load NREM epochs
load('D:\DPhil\Off-period detection\TestSignals\Data_Lukas\MH_20190217_DL_NA_VSspec.mat','nr')
%Find clean NREM segments (remove boundary epochs)
gaps=find(diff(nr)>1); %find last epoch of each episode
numepi=length(gaps); %number of NREM episodes (-1 because last epoch doesn't end with a gap)
cleanepochs=[];

%loop going through all NREM episodes (except for first and last)
for ep = 2:numepi
    
    Startepoch=nr(gaps(ep-1)+2); %find start epoch of this NREM episode
    Endepoch=nr(gaps(ep))-1; %find second last epoch of this NREM episode
    
    if Endepoch-Startepoch<3 %if NREM episode 2 epochs or shorter, there'll be no signal because w're cutting the first and last epoch
        continue
    else %if NREM episode at least 3 epochs, find bins corresponding to the start of the second epoch and end of the second last epoch
        cleanepochs=[cleanepochs; Startepoch Endepoch];
    end
end
numepochs=length(cleanepochs);

%Load LFP
load('D:\DPhil\Off-period detection\TestSignals\Data_Lukas\MH-LFP-20190217-ch8.mat','sig')
%Filter LFP
LFP_ch8=bandpass(sig,[0.5,30],fs_lfp);
clear sig

%load PNE
load('D:\DPhil\Off-period detection\TestSignals\Data_Lukas\pNe_data_MH_20190217_DL_Channels1to8','Ch8')
PNE_ch8=Ch8;
clear Ch8

%Extact NREM signals
%%% collect pNe signal for all NREM episodes
NremLFP=NaN(1,length(LFP_ch8)); %NaN vectors to be filled with NREM pNe signal
NremPNE=NaN(1,length(PNE_ch8));

%loop going through all NREM epochs (except for first and last)
for ep = 1:numepochs

Startepoch=cleanepochs(ep,1); %start of respective epoch
Endepoch=cleanepochs(ep,2); %end of respective epoch

%fill NaN vectors with LFP signal for this NREM episode
StartBinLFP=ceil(Startepoch*epochlength*fs_lfp);
EndBinLFP=floor(Endepoch*epochlength*fs_lfp);
NremLFP(StartBinLFP:EndBinLFP)=LFP_ch8(StartBinLFP:EndBinLFP);

%fill NaN vectors with pNe signal for this NREM episode
StartBinPNE=ceil(Startepoch*epochlength*fs_pne);
EndBinPNE=floor(Endepoch*epochlength*fs_pne);
NremPNE(StartBinPNE:EndBinPNE)=PNE_ch8(StartBinPNE:EndBinPNE);

end
%concatenate the signal from all selected NREM episodes to get rid of gaps
NremLFP=NremLFP(~isnan(NremLFP));
NremLFP=NremLFP*1000000;
NremPNE=NremPNE(~isnan(NremPNE));

% %Resample to pne fs for ease of plotting
% NremLFPresampled=resample(NremLFP,round(fs_pne),round(fs_lfp));

%%
%close all

LowDurThresh=42/1000;
MedDurThresh=102/1000;

secondsRange=floor(length(NremPNE)/fs_pne);
%startSec=randi(secondsRange);
%startSec=150;
endSec=startSec+8;
LFPsegment=NremLFP(fs_lfp*startSec:fs_lfp*endSec);
PNEsegment=NremPNE(fs_pne*startSec:fs_pne*endSec);
LFPtime=[startSec:1/fs_lfp:endSec]; 
PNEtime=[startSec:1/fs_pne:endSec];

startSegIN=find(OFF_periods_Lukas(8).OFFperiod(:,1)>fs_pne*startSec&OFF_periods_Lukas(8).OFFperiod(:,1)<fs_pne*endSec);
endSegIN=find(OFF_periods_Lukas(8).OFFperiod(:,2)>fs_pne*startSec&OFF_periods_Lukas(8).OFFperiod(:,2)<fs_pne*endSec);
SegIN=intersect(startSegIN,endSegIN);


figure('Position',[3.4,334.6,1526.4,290.4])
subplot(2,1,1)
plot(LFPtime,LFPsegment)
hold on
for i = 1:length(SegIN)
    startSegIN=(OFF_periods_Lukas(8).OFFperiod(SegIN(i),1)/fs_pne)*fs_lfp;
    endSegIN=(OFF_periods_Lukas(8).OFFperiod(SegIN(i),2)/fs_pne)*fs_lfp;
    plot(startSegIN*(1/fs_lfp):(1/fs_lfp):endSegIN*(1/fs_lfp),NremLFP(startSegIN:endSegIN),'r','LineWidth',2)     
    if i<length(SegIN)
        startSegIN2=(OFF_periods_Lukas(8).OFFperiod(SegIN(i+1),1)/fs_pne)*fs_lfp;
        if startSegIN2-endSegIN < LowDurThresh*fs_lfp
            plot(endSegIN*(1/fs_lfp):(1/fs_lfp):startSegIN2*(1/fs_lfp),NremLFP(endSegIN:startSegIN2),'k','LineWidth',1)     
        elseif LowDurThresh*fs_lfp<startSegIN2-endSegIN & startSegIN2-endSegIN < MedDurThresh*fs_lfp
            plot(endSegIN*(1/fs_lfp):(1/fs_lfp):startSegIN2*(1/fs_lfp),NremLFP(endSegIN:startSegIN2),'m','LineWidth',1)     
        end
    end
end
%ylim([-1000,1000])
subplot(2,1,2)
    plot(PNEtime,PNEsegment)
hold on
for i = 1:length(SegIN)
    startSegIN=OFF_periods_Lukas(8).OFFperiod(SegIN(i),1);
    endSegIN=OFF_periods_Lukas(8).OFFperiod(SegIN(i),2);
    plot(startSegIN*(1/fs_pne):(1/fs_pne):endSegIN*(1/fs_pne),NremPNE(startSegIN:endSegIN),'r')
end

%%
meanLFPampTotal=mean(NremLFP);
allLFPOFFsamples=[];
for i=1:length(OFF_periods_Lukas(8).OFFperiod)
    currentOFF=NremLFP((OFF_periods_Lukas(8).OFFperiod(i,1)/fs_pne)*fs_lfp:(OFF_periods_Lukas(8).OFFperiod(i,2)/fs_pne)*fs_lfp);
    allLFPOFFsamples=[allLFPOFFsamples,currentOFF];
end
meanLFPampOFF=mean(allLFPOFFsamples);

%%
%0.0205
%  402.9240


figure
subplot(2,1,1)
plot(LFPtime,LFPsegment)
hold on
for i = 1:length(SegIN)
    startSegIN=(OFF_periods_Lukas(8).OFFperiod(SegIN(i),1)/fs_pne)*fs_lfp;
    endSegIN=(OFF_periods_Lukas(8).OFFperiod(SegIN(i),2)/fs_pne)*fs_lfp;
    if endSegIN-startSegIN<=LowDurThresh*fs_lfp
        plot(startSegIN*(1/fs_lfp):(1/fs_lfp):endSegIN*(1/fs_lfp),NremLFP(startSegIN:endSegIN),'r','LineWidth',2)
    elseif LowDurThresh*fs_lfp<endSegIN-startSegIN & endSegIN-startSegIN<=MedDurThresh*fs_lfp
        plot(startSegIN*(1/fs_lfp):(1/fs_lfp):endSegIN*(1/fs_lfp),NremLFP(startSegIN:endSegIN),'m','LineWidth',2)
    else
        plot(startSegIN*(1/fs_lfp):(1/fs_lfp):endSegIN*(1/fs_lfp),NremLFP(startSegIN:endSegIN),'k','LineWidth',2)
    end   
end
subplot(2,1,2)
    plot(PNEtime,PNEsegment)
hold on
for i = 1:length(SegIN)
    startSegIN=OFF_periods_Lukas(8).OFFperiod(SegIN(i),1);
    endSegIN=OFF_periods_Lukas(8).OFFperiod(SegIN(i),2);
    if endSegIN-startSegIN<=LowDurThresh*fs_pne
        plot(startSegIN*(1/fs_pne):(1/fs_pne):endSegIN*(1/fs_pne),NremPNE(startSegIN:endSegIN),'r')
    elseif LowDurThresh*fs_pne<endSegIN-startSegIN & endSegIN-startSegIN<=MedDurThresh*fs_pne
        plot(startSegIN*(1/fs_pne):(1/fs_pne):endSegIN*(1/fs_pne),NremPNE(startSegIN:endSegIN),'m')
    else
        plot(startSegIN*(1/fs_pne):(1/fs_pne):endSegIN*(1/fs_pne),NremPNE(startSegIN:endSegIN),'k')
    end   
end    
