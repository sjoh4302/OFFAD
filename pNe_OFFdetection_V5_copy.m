%%%%
%{
GENERAL INFORMATION
This script compares different approaches to define an 'objective
threshold' in each individual pNe channel for subsequent analysis of
OFF-periods

AUTHOR INFORMATION
- script written by LBK and VVV in April 2020
    (contact:lukas.krone@dpag.ox.ac.uk)
- last update: LBK 20200515 (added pwelch approach but not 100% sure it's working correctly) 

INPUT
- mousename
- baseline date 'BL date'
- channel of interest 'chan'
- pNe signal 
- vigilance state info from scoring
- start and end time (in s) of time window (ideally 10 s) for 'test signal' 
with clear ON and OFF states during NREM sleep 

OUTPUT
- figure showing thresholds from differnt approaches on test signal 
%}
%% clear all variables, close all figures
clear all
close all

%% file and path description
path='D:\CorticalAreas\';
pathinPNE=[path,'OutputSignals24h\']; %path for pNe files
pathinVS=[path,'outputVSch\']; %path for VS files

mousename='MH';
BLdate='20190217';
chan='6';

SR=498.2462;

%provide time window of 10s NREM sleep as 'test signal' (unit seconds, find good trace in OpenScope)
TestWindowStart=9675;  %start time in sec; %TestWindowStart=round(TestWindowStart*SR); find respective pNe bin
TestWindowEnd=9685;  %end time in sec; %TestWindowEnd=round(TestWindowEnd*SR); find respective pNe bin
TestWindowLength=TestWindowEnd-TestWindowStart; %length in sec
% TestWindow=[TestWindowStart:TestWindowEnd];


%% load and process pNe signal
%{
1) convert it into absolute values
2) create resampled copy with new sampling rate of 256 Hz (for approach 2)
3) compare absolute signal before and after resampling
%}

%load pNe signal
filename=[mousename,'-pNe-',BLdate,'-ch',chan];
load([pathinPNE filename '.mat']);

%transform
sig=sig/1000000; %transform into uV values
sig=abs(sig); %take absolute values

%resample to 256 Hz
f1=SR;
f2=256;
num_chunks=4;
tail_length=20000;
visualize=0;

sig_resamp=resampling(sig,f1,f2,num_chunks,tail_length,visualize);


%plot figure to see resampling results
indSRorig=round(TestWindowStart*f1):1:round(TestWindowEnd*f1);
xOrig=10/length(indSRorig):10/length(indSRorig):TestWindowLength;
testsignal=sig(indSRorig);

indSRresamp=round(TestWindowStart*f2):1:round(TestWindowEnd*f2);
xResamp=10/length(indSRresamp):10/length(indSRresamp):TestWindowLength;
testsignal_resamp=sig_resamp(indSRresamp);

%check resampled signal visually
figure(1)
p1=plot(xOrig,testsignal,'b');
lgd_p1='original signal';
hold on
p2=plot(xResamp,testsignal_resamp,'r');
lgd_p2='resampled testsignal';

%figure settings
legend(lgd_p1,lgd_p2)
ylabel('Signal amplitude (in uV)')
xlabel('Time [s]')
title('test signal')
set(gcf,'Position',[200 200 1500 500])


%% collect pNe signal for all NREM episodes
% cutting off one 4 second epoch at the beginning and end of each episode to exclude transitional states

%load vigilance state information ('nr' are all artefact free NREM epochs)
filenameVS=[mousename,'-',BLdate,'-eeg17','-VSspec'];
load([pathinVS,filenameVS,'.mat'],'-mat','nr');

epochlength=4; %vigilance state scoring in 4s epochs

% find NREM episodes
gaps=find(diff(nr)>1); %find last epoch of each episode
numepochs=length(gaps); %number of NREM episodes (-1 because last epoch doesn't end with a gap)

%NaN vectors to be filled with NREM pNe signal
NremSig=NaN(1,length(sig));
NremSig_resamp=NaN(1,length(sig_resamp));

%loop going through all NREM epochs (exept for first and last)
for ep = 2:numepochs
    
    Startepoch=nr(gaps(ep-1)+1); %find start epoch of this NREM episode
    Endepoch=nr(gaps(ep))-1; %find second last epoch of this NREM episode
    
    if Endepoch-Startepoch<3 %if NREM episode 2 epochs or shorter, there'll be no signal because w're cutting the first and last epoch
        continue
    else %if NREM episode at least 3 epochs, find bins corresponding to the start of the second epoch and end of the second last epoch
        
        %fill NaN vectors with pNe signal for this NREM episode
        StartBin=ceil(Startepoch*epochlength*f1);
        EndBin=floor(Endepoch*epochlength*f1);
        NremSig(StartBin:EndBin)=sig(StartBin:EndBin);
        
        %fill NaN vectors with resampled pNe signal for this NREM episode
        StartBin_resamp=(Startepoch*epochlength*f2)+1; %first bin of second epoch
        EndBin_resamp=Endepoch*epochlength*f2; %last bin of second last epoch
        NremSig_resamp(StartBin_resamp:EndBin_resamp)=sig_resamp(StartBin_resamp:EndBin_resamp);
        
    end
end

totsig=NremSig(~isnan(NremSig));
totsig_resamp=NremSig_resamp(~isnan(NremSig_resamp)); %get rid of gaps


%% APPROACH 0 - visual determination of spike threshold
% based on test signal

figure(1)

%ask for visual threshold
prompt0 = 'optimal threshold based on visual inspection of test signal:';
threshold_0 = input(prompt0);

%plot visual threshold
hold on
p3=plot([0 TestWindowLength],[threshold_0 threshold_0],'k--','LineWidth',2);
lgd_p3='threshold based on inspection of test signal';
legend(lgd_p1,lgd_p2,lgd_p3)

%% APPROACH 1 - try to identify optimal cutoff in bimodal distribution of spike amplitudes or in smoothed test signal using different smoothing factors
%{
    1) plot histogram of spike amplitudes
    2) find cutoff point in bimodal distribution
    3) set threshold manually
%}

for av=1:16
    
    averagetotsig=movmean(totsig,av); %smooth total signal
    
    figure(10)
    subplot(4,4,av)
    h=histc(averagetotsig,[0:1:150]);
    bar(h);
    %     set(gca, 'XTick', [0:5:30]);
    %     xticklabels({[0:25:150]});
    xlabel('Spike Amplitude [uV]')
    
    averagetestsig=movmean(testsignal,av); %smooth test signal
    
    figure(11)
    subplot(4,4,av)
    plot(averagetestsig,'b')
end

%add figure titles
figure(10)
suptitle({['distribution of spike amplitutes in NREM signal'],['with increasing moving window size']})
figure(11)
suptitle({['smoothed test signal'],['with increasing moving window size']})

figure(1)

prompt1 = 'Threshold based on inspection of spike distribution (bimodality):';
threshold_1 = input(prompt1);

%plot visual threshold
hold on
p4=plot([0 TestWindowLength],[threshold_1 threshold_1],'g--','LineWidth',2);
lgd_p4='threshold based on inspection of spike distribution';
legend(lgd_p1,lgd_p2,lgd_p3,lgd_p4)


%% APPROACH 2 - determine threshold required the occurrence of 1 OFF period per second then move slightly up
% find threshold for 1 OFF period per second then move up by an arbitrary margin (e.g. 5% of the spike population)
%{
        1) find lowest absolute value in spike amplitude signal
        2) increase threshold stepwise by 1 to find threshold on which at least 1 OFF period (50 ms
        silence) per second is recorded
        3) add 5% more spikes from the population and use the amplitude of
        the highest spike as threshold
%}

%1 find lowest signal amplitude
[minVal,~]=min(totsig);
threshold_2=minVal; %start from lowest signal amplitude
GapsPerSec=0;

%2 stepwise threshold increase until OFF criterion fulfilled
TotsigLength=round(length(totsig)/SR); %length of total signal in sec


while GapsPerSec<=1
    threshold_2=threshold_2+1;
    
    NeuronalSpikes=find(totsig>threshold_2);
    Space=diff(NeuronalSpikes);
    LongGaps=find(Space > 50);
    
    GapsPerSec=length(LongGaps)/TotsigLength;
end

%move 5% up in the spike population
sortsig=sort(totsig);
threshold_percentile=find(sortsig==threshold_2);
threshold_percentile=threshold_percentile(end)/length(sortsig);

threshold_percentile_final=threshold_percentile+0.05; %add 5%
threshold_2_final=sortsig(round(length(sortsig)*threshold_percentile_final));

%plot signal sorted by amplitude
figure(20)
% semilogy(sortsig,'k','LineWidth',2)
plot(sortsig,'k','LineWidth',2)
hold on

%plot orignial threshold in distribution
cutoff=round(length(sortsig)*threshold_percentile);
plot([cutoff cutoff],[0 200]);
hold on

%plot final threshold in distribution
cutoff_final=round(length(sortsig)*threshold_percentile_final);
plot([cutoff_final cutoff_final],[0 200],'m:','LineWidth',2);
hold on

figure(1)
hold on
p5=plot([0 TestWindowLength],[threshold_2_final threshold_2_final],'m--','LineWidth',2);
lgd_p5='automatic threshold (1 OFF/sec + 5% more spikes)';
legend(lgd_p1,lgd_p2,lgd_p3,lgd_p4,lgd_p5)

%% APPROACH 3 - find peak of the power spectral density (PSD) estimate of the NREM pNe signal, find threshold so that the frequency of OFF periods matches the peak frequency of the slow osciallation 

%{
1) create PSD estimate for the pNe signal
2) increase threshold stepwise by 1 to find threshold on which the frequency of OFF periods (>50 ms silence) matches the frequency of the slow oscillation

%}

%%% use pwelch to create PSD (CAVE: I'm not sure this is correct)
signal = sig_resamp;
window = 4*f2; %4 seconds * sampling rate of signal
noverlap = window/2; %half window overlap recommended in documentation
f=0:0.25:30; %my frequencies of interest 
fs=f2; %sampling rate 

[pxx,f] = pwelch(signal,window,noverlap,f,fs);

figure(30)
plot(f,10*log10(pxx))
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')

%find maximum of pWelch >0.5 Hz 
[~,bin]=max(pxx(3:end))
PeakFreq=f(2+bin);
NumOFFperSec=PeakFreq;


%1 find lowest signal amplitude
[minVal,~]=min(totsig);
threshold_3=minVal; %start from lowest signal amplitude
GapsPerSec=0;

%2 stepwise threshold increase until OFF criterion fulfilled
TotsigLength=round(length(totsig)/SR); %length of total signal in sec

while GapsPerSec<=NumOFFperSec
    threshold_3=threshold_3+1;
    
    NeuronalSpikes=find(totsig>threshold_3);
    Space=diff(NeuronalSpikes);
    LongGaps=find(Space > 50);
    
    GapsPerSec=length(LongGaps)/TotsigLength;
end

figure(1)
hold on
p6=plot([0 TestWindowLength],[threshold_3 threshold_3],'c--','LineWidth',2);
lgd_p6='automatic threshold (derived from spectral peak of pNe signal)';
legend(lgd_p1,lgd_p2,lgd_p3,lgd_p4,lgd_p5,lgd_p6)



%% APPROACH 4 - optimise number of OFF periods

OFFnum=NaN(1,200);
thresh=0;

for thr = 1:200
    
    thresh=thresh+1;
    NeuronalActivity=find(totsig>thresh);
    Space=diff(NeuronalActivity);
    LongGaps=find(Space > 50);
    
    OFFnum(thr)=length(LongGaps);
end

figure(40)
plot(OFFnum)
title('Number of OFF periods for all thresholds')

[y,threshold_4]=max(OFFnum); %lowest threshold on which OFF periods are max

figure(1)
p7=plot([0 TestWindowLength],[threshold_4 threshold_4],'y--','LineWidth',2);
lgd_p7='automatic threshold (based on maximum number of OFF-periods)';
legend(lgd_p1,lgd_p2,lgd_p3,lgd_p4,lgd_p5,lgd_p6,lgd_p7)




%% OFF period detection

%%%% continue here! 
%{
Step 1: ask Vlad which method shall be used for thresholding
Step 2: discuss how to define OFF periods based on the threshold 
  Vlad's old approach: take all 'gaps' between 50 and 4000 ms, take duration for top
  (longest) 20% on baseline day as threshold and compare duration and
  number in respect to sleep homeostasis 
  other ideas: 
  - define ON periods (burst firing on either side (maybe 5 spikes in 50
  ms or so) and define gaps < 4000ms  between them as OFF periods, even if individual spikes occur
  (spikes during OFF period can be used for 'OFF period quality' determination)
  - only select gaps with positive LFP waves
%}

