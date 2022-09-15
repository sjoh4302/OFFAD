function [OFFDATA]=OFFAD_clustering(OFFDATA)
%
% Clustering page: Perform clustering for OFF period detection 
%
% Author: Christian Harding 2022
% OFF Period Automated Detection (OFFAD) toolbox
% christian.harding@sjc.ox.uk
%
% Requires:
% - Statistics and Machine learning toolbox
% - Signal processing toolbox
%


% Start timer to compute clustering length
functionTimer = tic;

% Get screensize to set figure position
screensize = get( groot, 'Screensize' );   
 
% Create figure for clustering warning
g.Clustering = figure('name', 'OFFAD (OFF_period Automated Detection) - Clustering', ... 
	'numbertitle', 'off', ...
	'Position',[screensize(3)/6,3*screensize(4)/8,4*screensize(3)/6,screensize(4)/4], ...
    'Toolbar','none',...
    'Menubar','none',...
	'Tag','OFFAD_CLUSTER');
uicontrol(g.Clustering,'Style', 'text','String','Warning: Clustering in progress',...
    'FontWeight','bold','FontSize',22,...
    'Units','normalized',...
    'Position',[0.2 0.2 0.6 0.6]);

drawnow 

%%%%% Perform clustering

% Step 1: Extract NREM epochs to train GMM

% Extract NREM episodes (consecutive NREM epochs)
vigState = load(OFFDATA.VSpathin,'-mat','nr');
vigState = vigState.('nr');

% Remove epochs at the beginning and end of each episode to exclude transitional states
endEpi=find(diff(vigState)>1); %find last epoch of each episode
endEpi=[endEpi;length(vigState)]; %add final episode end
startEpi=find(diff(vigState)>1)+1; %find first epoch of each episode
startEpi=[1;startEpi]; %add first episode start
numepi=length(startEpi); %number of NREM episodes
cleanepochs=[];

% Loop going through all NREM episodes to find start and end epoch times
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


% Intialise structure to store clustering results
exampleObject = matfile(OFFDATA.MUApathin);
MUAlength=size(exampleObject,OFFDATA.ChannelsFullName(1),2);
OFFDATA.StartOP=sparse(repmat(logical(0),MUAlength,length(OFFDATA.ChannelsFullName)));
OFFDATA.EndOP=sparse(repmat(logical(0),MUAlength,length(OFFDATA.ChannelsFullName)));
OFFDATA.AllOP=sparse(repmat(logical(0),MUAlength,length(OFFDATA.ChannelsFullName)));
clear exampleObject 


% Repeat extraction for each channel
for chanNum = 1:length(OFFDATA.ChannelsFullName)

    % Step 2: Extract MUA signal from NREM epochs for training
    %{
    1) Convert it into absolute values
    2) Extract and concatenate NREM episodes only
    %}


    % Start channel clustering timer
    channelTimer=tic;


    chan=OFFDATA.ChannelsFullName(chanNum);    
    display(['Clustering ',char(chan)])

    % Load MUA signal for selected channel
    MUA = load(OFFDATA.MUApathin,'-mat',chan);
    MUA = MUA.(chan);


    % Convert to uV
    if OFFDATA.MUAunit=="uV*1000000" 
        MUA=MUA/1000000;
    end

    % Convert to lower memory data class if possible
    if class(MUA)=="single"
        MUA = int16(MUA);
    end

    % Return warning if incorrect amplitude detected
    if length(unique(MUA))<4
        warning('Wrong amplitude units specified')
    end

    % Convert to absolute values
    MUA = abs(MUA);

    % Create empty variables to store MUA signals from NREM only
    vsMUA=NaN(1,length(MUA)); %NaN vectors to be filled with MUA signal
    vsMUAtime=NaN(1,length(MUA)); %NaN vectors to be filled with recording time values

    % Loop going through all NREM epochs (except for first and last)and extract MUA
    for ep = 1:numepochs

        Startepoch=cleanepochs(ep,1); %start of respective epoch
        Endepoch=cleanepochs(ep,2); %end of respective epoch

        %fill NaN vectors with MUA signal for this NREM episode
        StartBin=ceil(Startepoch*OFFDATA.epochLen*OFFDATA.MUAfs);
        EndBin=floor((Endepoch-1)*OFFDATA.epochLen*OFFDATA.MUAfs);
        try
            vsMUA(StartBin:EndBin)=MUA(StartBin:EndBin);
            vsMUAtime(StartBin:EndBin)=StartBin:EndBin;
        end

    end
    % Concatenate the signal from all selected NREM episodes to get rid of gaps
    vsMUA=vsMUA(~isnan(vsMUA));
    vsMUAtime=vsMUAtime(~isnan(vsMUAtime));


    % Step 3: Extract MUA signal from all vigialnce states for clustering
    
    % Create variables to store MUA
    cleanMUA=NaN(1,MUAlength); %NaN vectors to be filled with MUA signal
    cleanMUAtime=NaN(1,MUAlength); %NaN vectors to be filled with recording time values

    % Loop going through all epochs from all states (except for first and last)and extract MUA
    states=["nr","w","r","mt"];
    for st=1:length(states)
        vigState = load(OFFDATA.VSpathin,'-mat',states(st));
        vigState = vigState.(states(st));

        for k=1:length(vigState)
            startVigEpoch=floor((vigState(k)-1)*OFFDATA.epochLen*OFFDATA.MUAfs)+1;
            endVigEpoch=floor((vigState(k))*OFFDATA.epochLen*OFFDATA.MUAfs);

            if endVigEpoch>MUAlength
                warning('Epoch skipped. MUA shorter than expected. Check sampling rate is correct')
            else
                cleanMUA(startVigEpoch:endVigEpoch)=MUA(startVigEpoch:endVigEpoch);
                cleanMUAtime(startVigEpoch:endVigEpoch)=1;
            end             
        end
    end
    % Concatenate the signal from all selected episodes to get rid of gaps
    cleanMUA=cleanMUA(~isnan(cleanMUA));
    cleanMUAtime=find(cleanMUAtime==1);



    % Step 4: Apply smoothing windows to MUA signal

    % Apply heavy smoothing (Guassian window) to MUA signals
    if OFFDATA.clustVar1Select==1
        winSize=(OFFDATA.clustVar1Smooth*2)+1;
        w1=gausswin(winSize);
        clusterVar1=conv(vsMUA,gausswin(winSize))/sum(w1);
        clusterVar1_full=conv(cleanMUA,gausswin(winSize))/sum(w1);
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
        [clusterVar1,F] = offPeriodScalogram(vsMUA,OFFDATA.MUAfs,percOverlap,LB_freq);
        [clusterVar1_full,F] = offPeriodScalogram(cleanMUA,OFFDATA.MUAfs,percOverlap,LB_freq);
        clear LB_freq percOverlap F
    end

    % Apply light smoothing (Guassian window) to MUA signals
    if OFFDATA.clustVar2Select==1
        winSize=(OFFDATA.clustVar2Smooth*2)+1;
        w2=gausswin(winSize);
        clusterVar2=conv(vsMUA,gausswin(winSize))/sum(w2);
        clusterVar2_full=conv(cleanMUA,gausswin(winSize))/sum(w2);
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
        [clusterVar2,F] = offPeriodScalogram(vsMUA,OFFDATA.MUAfs,percOverlap,LB_freq);
        [clusterVar2_full,F] = offPeriodScalogram(cleanMUA,OFFDATA.MUAfs,percOverlap,LB_freq);
        clear LB_freq percOverlap F
    end

    clear MUA
    clear vsMUA

    
    
    % Step 5: Fit Gaussian mixture model to establish components
    
    % Create try/catch loop in case of ill-conditioning errors
    try
        % If optimal cluster number not established already in preclustering phase, run
        % preclustering  
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
        
        % Collect smoothed data
        allData=[clusterVar1_full',clusterVar2_full'];
        % Select subsample of data to fit GMM 
        randPoints=randi(length(clusterVar1),round(length(clusterVar1)*(OFFDATA.percSamp/100)),1);
        sampData=[clusterVar1(randPoints)',clusterVar2(randPoints)'];
        options = statset('MaxIter',1000,'TolFun',1e-5);

        % Fit GMM to smoothed data (if models not specified in Preset mode)
        if isempty(OFFDATA.GMModels.(OFFDATA.ChannelsFullName(chanNum)))
            GMModel = fitgmdist(sampData,OFFDATA.OptimalK(chanNum),'Replicates',20,'Options',options);
            OFFDATA.GMModels.(OFFDATA.ChannelsFullName(chanNum)) = GMModel; 
        else
            GMModel = OFFDATA.GMModels.(OFFDATA.ChannelsFullName(chanNum));
        end
        
        % Step 6: Cluster MUA using posterior probabilities of GMM components
        
        % Cluster all data
        IDX = cluster(GMModel,allData);
        
        % Identify points beloinging to the smallest amplitude cluster 
        % (i.e. the OFF period cluster)
        [~,offCluster]=min(GMModel.mu(:,1));
        OFF_clust_points=cleanMUAtime(IDX==offCluster);

        % Step 7: Identify OFF periods as MUA negative half waves containing OFF period points 
        
        % Create empty variable to store raw OFF period points
        clustOP=sparse(OFF_clust_points,1,logical(1),length(OFFDATA.StartOP),1,length(OFF_clust_points));
        
        % Load MUA signal
        MUA = load(OFFDATA.MUApathin,'-mat',chan);
        MUA = MUA.(chan);

        % Convert to uV
        if OFFDATA.MUAunit=="uV*1000000" 
            MUA=MUA/1000000;
        end

        % Produce warning is wrong amplitude suspected
        if length(unique(MUA))<4
            warning('Wrong amplitude units specified')
        end
        
        % Convert to absolute values
        MUA = abs(MUA); %take absolute values 
        
        % Find baseline amplitude (mean MUA amplitude during WAKE)
        if isnan(OFFDATA.BaselineAmp(chanNum))
            
            % Load vigilance state information for this animal and recording date ('nr' are all artefact free NREM epochs)
            wakeState = load(OFFDATA.VSpathin,'-mat','w');
            wakeState = wakeState.('w');

            % Find WAKE episodes
            endEpi=find(diff(wakeState)>1); %find last epoch of each episode
            endEpi=[endEpi;length(wakeState)]; %add final episode end
            startEpi=find(diff(wakeState)>1)+1; %find first epoch of each episode
            startEpi=[1;startEpi]; %add first episode start
            numepi=length(startEpi); %number of NREM episodes
            wakeEpochs=[];

            %Find all WAKE epochs 
            for ep = 1:numepi

                Startepoch=wakeState(startEpi(ep)); %find start epoch of this NREM episode
                Endepoch=wakeState(endEpi(ep)); %find second last epoch of this NREM episode

                if Endepoch-Startepoch<2 %if NREM episode 2 epochs or shorter, there'll be no signal because we're cutting the first and last epoch (i.e. state transitions)
                    continue
                else %if NREM episode at least 3 epochs, find bins corresponding to the start of the second epoch and end of the second last epoch
                    wakeEpochs=[wakeEpochs; Startepoch Endepoch];
                end
            end

            wakeMUA=NaN(1,length(MUA)); %NaN vectors to be filled with MUA signal

            % Loop going through all wake epochs (except for first and last)
            for ep = 1:length(wakeEpochs)

                Startepoch=wakeEpochs(ep,1); %start of respective epoch
                Endepoch=wakeEpochs(ep,2); %end of respective epoch

                %fill NaN vectors with MUA signal for this NREM episode
                StartBin=ceil(Startepoch*OFFDATA.epochLen*OFFDATA.MUAfs);
                EndBin=floor((Endepoch-1)*OFFDATA.epochLen*OFFDATA.MUAfs);
                try
                    wakeMUA(StartBin:EndBin)=MUA(StartBin:EndBin);
                end
            end
            
            % Concatenate the signal from all selected WAKE episodes to get rid of gaps
            wakeMUA=wakeMUA(~isnan(wakeMUA));
            % Compute WAKE mean
            wakeMUAmean=mean(wakeMUA);
            OFFDATA.BaselineAmp(chanNum)=wakeMUAmean;
            clear wakeMUA wakeEpochs
        else
            wakeMUAmean = OFFDATA.BaselineAmp(chanNum);
        end

        % Subtract baseline amplitude from MUA
        relativeMUA=MUA-wakeMUAmean;
       
        % Find MUA below baseline 
        polarityMUA=ones(length(relativeMUA),1);
        polarityMUA(relativeMUA<0)=-1;
        
        % Extract MUA between negative zero-crossings
        nhwStarts=find(diff(polarityMUA)<0)+1; %NHW=negative half waves
        nhwEnds=find(diff(polarityMUA)>0);
        if nhwEnds(1)<nhwStarts(1)
            nhwEnds=nhwEnds(2:end);
        end
        if nhwStarts(end)>nhwEnds(end)
            nhwStarts=nhwStarts(1:end-1);
        end
        
        % Find intersect between negative zero-crossings and clustered OFF
        % period points
        OFFperiodIndex=logical([]);
        for i = 1:length(nhwStarts)
            if sum(clustOP(nhwStarts(i):nhwEnds(i)))>0
                OFFperiodIndex(i)=logical(1);
            else
                OFFperiodIndex(i)=logical(0);
            end
        end
        clear MUA
   
        % Start times of negative zero-crossings containing OFF period points
        StartOFF=nhwStarts(OFFperiodIndex);
        % End times of negative zero-crossings containing OFF period points
        EndOFF=nhwEnds(OFFperiodIndex);
        % All OFF period times
        tempallOP=zeros(length(OFFDATA.StartOP),1);
        for i=1:length(StartOFF)
            tempallOP(StartOFF(i):EndOFF(i))=1;
        end

        AllOFF=find(tempallOP==1);
        clear tempBorders newCounter tempallOP

    catch
        warning('Failed to cluster channel')
        StartOFF=[];
        EndOFF=[];
        AllOFF=[];

    end

    % Store START OFF-P data
    OFFDATA.StartOP(:,chanNum)=sparse(StartOFF,1,logical(1),length(OFFDATA.StartOP),1,length(StartOFF));

    % Store END OFF-P data
    OFFDATA.EndOP(:,chanNum)=sparse(EndOFF,1,logical(1),length(OFFDATA.StartOP),1,length(EndOFF));

    % Store ALL OFF-P data
    OFFDATA.AllOP(:,chanNum)=sparse(AllOFF,1,logical(1),length(OFFDATA.StartOP),1,length(AllOFF));


    clear OFFperiod OFF_clust_points ON_clust_points clusterVar1 clusterVar2 ...
          clusterIDX_total clustVar1ThreshScaled clustVar2ThreshScaled StartOFF EndOFF ...
          nhwStarts nhwEnds OFFperiodIndex AllOFF

    % Display total channel clustering time
    channelTime=toc(channelTimer);
    display(['Channel clustering time = ',char(string(channelTime)),' sec'])
end

% Display total function clustering time
functionTime=toc(functionTimer);
display(['Total clustering time = ',char(string(functionTime)),' sec'])

% Return to main menu when page deleted    
delete(findobj('Tag','OFFAD_IMPORT'));
OFFAD('redraw');
set(findobj('Tag','OFFAD_CLUSTER'),'Visible','Off');

