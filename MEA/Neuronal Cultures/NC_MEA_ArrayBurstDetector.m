function ArrayBurstData = NC_MEA_ArrayBurstDetector(MEAData)
%This function is meant to be used to detect Bursts in recordings using the
%multichannel systems 60 channel MEA. 
%
%Inputs
%   MEAData           - Data
%   ChannelNumber     - Desired Channel Number to get the Burst Parameters
%
%Outputs
%   BurstData          - Data about the Bursts
%       BurstLength    = Average burst length
%       SpikesPerBurst = Average amount of spikes per burst
%       BurstFreq      = Number of Burst that occur in the recording 
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN:  April 12, 2018
%LAST MODIFIED: April 18, 2018
%Modified from MEA_BurstParamCalc
%v1.0

%Define the desired Thresholds
BurstSpikeThreshold = 3;  %Number of Spikes per Burst
MinNumChannels      = ActiveChannels*0.10; %Minumun number of channels that
                                           %checked in order for a burst to
                                           %be considered an array burst                                    %burst. 
if MEAData.PeakData(ChannelNumber).SpikeNum > 2
    %Load the desired data
    SpikeTimes = MEAData.PeakData(ChannelNumber).SpikePeakTime;
    ISI        = MEAData.PeakData(ChannelNumber).ISI;
    BurstNum   = 0;
    counter  = 1;
    for ii = 1:(length(ISI)-1)
        
        
        
        %Calculate the difference in the interspike intervals
        DeltaISI = ISI(ii+1) - ISI(ii);
        
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
                Burst(BurstNum).BurstLength    =  SpikeTimes(ii) - SpikeTimes((ii - counter) + 1);
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
        BurstData.FrequencyofBurst = 0;
        BurstData.AveBurstLength   = 0;
        BurstData.AveSPB           = 0;
        BurstData.NumBurst         = 0;
    %If there is at least one burst, define all burst parameters
    elseif BurstNum >= 1
        %determine the length of the recording
        RecordingLength = MEAData.Time(ChannelNumber).RecordingTime;
        
        
        for ii = 1:BurstNum
            BurstData.Burst(ii).BurstLength      = Burst(ii).BurstLength;
            BurstData.Burst(ii).SpikesPerBurst   = Burst(ii).SpikesperBurst;
            %Make Burst parameters into vector to calculate the averages
            TempBL(ii)  = Burst(ii).BurstLength;
            TempSPB(ii) = Burst(ii).SpikesperBurst;
        end
        
        
        BurstData.FrequencyofBurst = BurstNum / RecordingLength;
        BurstData.AveBurstLength = mean(TempBL);
        BurstData.AveSPB         = mean(TempSPB);
        BurstData.NumBurst       = BurstNum;
    else
        fprintf('ERROR\n **REVISE CODE1**')
    end
else
    fprintf('ERROR\n **REVISE CODE Chan %d**', ChannelNumber)
    
end