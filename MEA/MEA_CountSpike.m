function MEA_CountSpike(MEAData, NumFiles, NumChannels, ShowDetectedPeaks)
%This is meant to show how the spike rate changes overtime in each channel.
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: July 18, 2017
%LAST MODIFIED: July 24, 2017
%v1.0

%Define any input not specified to their default value
if ~exist('ShowDetectedPeaks','var')
     %if ShowDetectedPeaks does not exist, default set it to N
      ShowDetectedPeaks = 'N';
      
end

%This part is where the number of spikes for each recording is determined.
%This section is the most time consuming and is dependent on the number of
%channels and files. This might be a good place to find ways to be more
%efficient
NumSpikesFile = struct;
tic 
for ii = 1:NumFiles
    TempNumSpikes = M_SpikeCounter(MEAData(ii), NumChannels, ShowDetectedPeaks);
    for n = 1:length(TempNumSpikes)
        NumSpikesFile(ii).NumSpikesChannel(n) = TempNumSpikes(n);
    end
end
toc
%Creates a structure used in the rest of the code to easily and logically
%access the data of each channel thorugh time (throught the different
%recordings)
Channel = struct;
for ii = 1:NumChannels
    for n = 1:NumFiles
    Channel(ii).File(n) =  NumSpikesFile(n).NumSpikesChannel(ii);
    end
end

%Determine which channels not to plot
%set a minimum number of spikes each 2 min recording has to have in order
%to be considered active and set the number of recordings that need to be 
%active in order to classify a channel to be active throughout time
P_SpikeThresh = 20;
P_FileThresh  = NumFiles*0.25;

CPurged = [];
Eliminated = 0;
%Calculate which Channels have less than Spike threshold 
for ii = 1:NumChannels
    PurgeVal = 0;
    for iii = 1:NumFiles
        if Channel(ii).File(n) < P_SpikeThresh
            PurgeVal = PurgeVal + 1;
        end
    end
    %eliminate the channel if it did not meet the criteria
    if PurgeVal > P_FileThresh
        Eliminated = Eliminated + 1;
        CPurged(Eliminated) = ii;
    end
    
end

%Print the channels that were eliminated
fprintf('Eliminated the following channels due to low activity:\n')
for ii = 1:length(CPurged)
    Num = MEAData(1).Settings.InfoChannel.Label{CPurged(ii)};
    fprintf('Channel %s \n', Num);

end
fprintf('Eliminated %d channels due to low activity\n', Eliminated)

%Determine which channels did not experience any change in firing frequency
EliminatedAct = 0;

%Set the Threshold for the channels to have experineced change in their
%firing rate compared tot he first recording (Baseline)
ChangeThreshUpper = 1.1;
ChangeTreshLower  = 0.90;
FileChangeThresh = NumFiles*0.80;

fprintf('Eliminated the following channels due small change in activity\n')
for ii = 1:NumChannels
    %do not look at channels already eliminated due to low spiking activity
    if ismember(ii, CPurged) == 0
        PurgeVal = 0;
        for iii = 2:NumFiles
            FoldChange = (Channel(ii).File(iii))*(1/(Channel(ii).File(1)));
            %determine if the channel had a decrease or increase in
            %activity throught the recordings
            if  FoldChange > ChangeThreshUpper || FoldChange < ChangeTreshLower
            else
                PurgeVal = PurgeVal + 1;
            end
        end
        %Eiminate the channel if it did not meet the criteria
        if PurgeVal > FileChangeThresh 
            Eliminated = Eliminated + 1;
            EliminatedAct = EliminatedAct + 1;
            CPurged(Eliminated) = ii;  
            Num = MEAData(1).Settings.InfoChannel.Label{ii};
            fprintf('Channel %s \n', Num)
        end
    end
end

fprintf('Eliminated %d channels due to no change in activity levels', EliminatedAct)

%Make the Total Spike plot
figure
hold on  
%Change displayed names
ChannelNames = cell(1,NumChannels-length(CPurged));
plotted = 0;  

for ii = 1:NumChannels

    if ismember(ii, CPurged) == 0
    %Define the axis values
    x = 1:NumFiles;
    y = Channel(ii).File;
    plot(x,y)

    ylabel('Number of Spikes')
    title('Spikes per channel')

    plotted = plotted+1;
    
    ChannelLabel = strcat('Channel-', MEAData(1).Settings.InfoChannel.Label{ii});
    ChannelNames(plotted) = {ChannelLabel};
    
    end
end

legend(ChannelNames)
legend('hide')
hold off

figure
hold on
% For changing display names
ChannelNames = cell(1,NumChannels-length(CPurged));
%Normalized Spike plot
plotted = 0; 

for ii = 1:NumChannels
    
    if ismember(ii, CPurged) == 0
    
    %Define the axis values
    x = 1:NumFiles;
    
    y = (Channel(ii).File)*(1/(Channel(ii).File(1)));
    plot(x,y)
    
    ylabel('Percent Baseline')
    title('Normalized Spikes per channel')
    
    %Change displayed names
    
    plotted = plotted+1;
    
    ChannelLabel = strcat('Channel-', MEAData(1).Settings.InfoChannel.Label{ii});
    ChannelNames(plotted) = {ChannelLabel};
    
    end
end

legend(ChannelNames)
legend('hide')
end

    