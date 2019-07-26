function NumSpikes = M_SpikeCounter(MEAData, NumChannels, ShowDetectedPeaks)
%This function is used to determine how many spikes occured in each channel
%of a MEA recording
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: July 18, 2017
%LAST MODIFIED: July 18, 2017
%v1.0


%
Channels = 1:NumChannels; 
%Preallocate for speed and memory
SpikesD = struct;


%Get the Spike data for each channel and save it to the struct
for ii = 1:NumChannels
   Temp_PeakData = MEA_Spike_Detector(MEAData, Channels(ii), ShowDetectedPeaks);
   %Define hte parameters in the struct
   SpikesD.PeakData(ii)    = Temp_PeakData;
      
end
NumSpikes = zeros(1, NumChannels);
for ii = 1:NumChannels
    NumSpikes(ii) = length(SpikesD.PeakData(ii).SpikePeakTime);
end