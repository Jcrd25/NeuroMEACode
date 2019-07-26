function ArrayBurstData = NC_ArrayBurstDetector(MEANCData)
%This function detects Array or netwrok bursting activity. It is intended
%to be used for 60 microelectrode arrays. 
%Input
%   MEAData - has the electrophysiological data from all the channels of
%   the microeletreode array
%Output
%   FileData - Summarized version of all the time recordings of a single
%   slice
%Outputs
%   BurstData          - Data about the Bursts
%       BurstLength    = Average burst length
%       SpikesPerBurst = Average amount of spikes per burst
%       BurstFreq      = Number of Burst that occur in the recording 
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: February 11, 2019
%LAST MODIFIED: 
%v1.0

NumChannels = length(MEANCData.PeakData);
%Define the desired Thresholds
ISIThreshold        = 0.150; %(in s)
BurstSpikeThreshold = 12;  %Number of Spikes per Burst
SpikeFreqThreshold  = 0.2*MEANCData.Time(1).RecordingTime;

%predefine the active channel variables
ActiveChanNum = 0;
ActiveChan = zeros(NumChannels,1);

%Active electrodes show firing rates of at least 0.2 Hz
for ChanNum = 1:NumChannels
    if MEANCData.PeakData(ChanNum).SpikeNum >= SpikeFreqThreshold
        ActiveChan(ChanNum) = true;
        ActiveChanNum = ActiveChanNum+1;
    else
        ActiveChan(ChanNum) = false;
    end
end
%A minimum of 4 channels need to be actively participating in the burst for
%it to be considered an array burst
if ActiveChanNum < 20
    ChanBurstThreshold  = 4;
else
    %Need at least 20% of the channels to be active in the bursting
    ChanBurstThreshold  = ActiveChanNum*.20; 
end
%Concatenate all of the spikepeaktimes into one vector.
SpikePeakTimesConc =[];
SpikePeakLocConc =[];
for ii = 1:NumChannels
    if ActiveChan(ii) == 1
        SpikePeakTimesConc = [SpikePeakTimesConc MEANCData.PeakData(ii).SpikePeakTime];
        SpikePeakLocConc   = [SpikePeakLocConc MEANCData.PeakData(ii).SpikePeakLoc'];
    end
end
SpikePeakTimesConc = sort(SpikePeakTimesConc);
SpikePeakLocConc   = sort(SpikePeakLocConc);



%Check if the Channel can be analyzed for burst data.
if length(SpikePeakTimesConc) <= 2
    %if Burst cannot be calculated then define all parameters as zero
    ArrayBurstData.FrequencyofBurst = 0;
    ArrayBurstData.AveBurstLength   = 0;
    ArrayBurstData.AveSPB           = 0;
    ArrayBurstData.NumBurst         = 0;
    
%Due to the definition of a burst used, we need to have at least more than
%2 spikes
elseif length(SpikePeakTimesConc) > 2
    BurstNum   = 0;
    counter  = 1;
    for ii = 1:(length(SpikePeakTimesConc)-1)
        
        
        
        %Calculate the difference in the interspike intervals
        DeltaISI = SpikePeakTimesConc(ii+1) - SpikePeakTimesConc(ii);
        
        %Determine if the ISI are within the desired closeness
        if DeltaISI <= ISIThreshold
            counter = counter + 1;
        %Once there are no longer consecutive close enough spikes,
        %determine if there are enough spikes for it to be consider a
        %Burst
        elseif DeltaISI > ISIThreshold
            %Define the Burst properties if it is a burst 
            if counter >= BurstSpikeThreshold
                BurstNum = BurstNum+1;
                Burst(BurstNum).Onset          = SpikePeakTimesConc((ii - counter) + 1);             
                Burst(BurstNum).End            = SpikePeakTimesConc(ii);
                Burst(BurstNum).OnsetLoc       = SpikePeakLocConc((ii - counter) + 1);             
                Burst(BurstNum).EndLoc         = SpikePeakLocConc(ii);
                Burst(BurstNum).BurstLength    = SpikePeakTimesConc(ii) - SpikePeakTimesConc((ii - counter) + 1);
                Burst(BurstNum).SpikesperBurst = counter+1;
                %reset counter
                counter = 1;
            %If there are not enough spikes to clasify it as a
            %Burst reset the counter
            elseif counter < BurstSpikeThreshold
                counter = 1;
            end
        end
    end
    
    %Define BurstData
    
    %If there are no burst, all burst parameters are considered 0
    if BurstNum == 0
        ArrayBurstData.FrequencyofBurst = 0;
        ArrayBurstData.AveBurstLength   = 0;
        ArrayBurstData.AveSPB           = 0;
        ArrayBurstData.NumBurst         = 0;
    %If there is at least one burst, define all burst parameters
    elseif BurstNum >= 1
        %determine the length of the recording
        RecordingLength = MEANCData.Time(1).RecordingTime;
        
        
        for ii = 1:BurstNum
            ArrayBurstData.Burst(ii).Onset            = Burst(ii).Onset;
            ArrayBurstData.Burst(ii).End              = Burst(ii).End;
            ArrayBurstData.Burst(ii).OnsetLoc         = Burst(ii).OnsetLoc;
            ArrayBurstData.Burst(ii).EndLoc           = Burst(ii).EndLoc;
            ArrayBurstData.Burst(ii).BurstLength      = Burst(ii).BurstLength;
            ArrayBurstData.Burst(ii).SpikesPerBurst   = Burst(ii).SpikesperBurst;
            %Make Burst parameters into vector to calculate the averages
            TempBL(ii)  = Burst(ii).BurstLength;
            TempSPB(ii) = Burst(ii).SpikesperBurst;
        end
        
        
        ArrayBurstData.FrequencyofBurst = BurstNum / RecordingLength;
        ArrayBurstData.AveBurstLength = mean(TempBL);
        ArrayBurstData.AveSPB         = mean(TempSPB);
        ArrayBurstData.NumBurst       = BurstNum;
    else
        fprintf('ERROR\n **REVISE CODE1**')
    end
else
    fprintf('ERROR\n **REVISE CODE Chan %d**', ChannelNumber)
    
end

