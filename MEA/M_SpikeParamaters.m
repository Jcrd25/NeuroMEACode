function M_SpikeParamaters(MEAData)
%This function is used to extract the spikes detetected and save them in a
%struct.
%
%Input
%Output
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: September 20, 2017
%LAST MODIFIED: September 20, 2017
%v1.0

%Detect the spikes and get the peak time stamp
%Get the Spike data for each channel and save it to the struct
for ii = 1:NumChannels
   Temp_PeakData = MEA_Spike_Detector(MEAData, Channels(ii), ShowDetectedPeaks);
   %Define hte parameters in the struct
   SpikesD.PeakData(ii)    = Temp_PeakData;
      
end
%Define all the variables required
%Define the time before and after the peak that will be considered to be a
%spike
SpikeBeg = 3; %in ms
SpikeEnd = 3; %in ms

SpikeNum = length(SpikeD.PeakData)
%Frequency of recording in Hz
Fs = 1/SamplingRate;
%How much time the recording takes
TimeEnd = length(MEAData)*Fs;
%The time values of the data 
Time = 0:Fs:TimeEnd;
PreSpikeTime  = SpikeBeg/Fs;
PostSpikeTime = SpikeEnd/Fs; 
%Count how many spikes where extracted and which were not
    SpikesExtracted    = 0;
    SpikesNotExtracted = 0;
   
    %For 1st spike   
    
   for n = 1:NumSpikes
       %Extract Spikes using the PreSpikeTime. The spike is extracted from
       %the amount of time specified by PreSpikeTime before the Peak until the next spike starts.
       %It for the last spike see below.
       if n == 1
           %determine if the first spike peak has at least 330 points before, if not
           %do not use it
           DeltaTime1Spike = PeakLocation(n)+1;
           if DeltaTime1Spike > PreSpikeTime
               SpikesData.Spike(n).SpikeTime = Time(PeakLocation(n)-PreSpikeTime : PeakLocation(n) + PostSpikeTime);
               SpikesData.Spike(n).SpikeVoltage = Voltage(PeakLocation(n)-PreSpikeTime : PeakLocation(n+1) + PostSpikeTime);
               SpikesData.Spike(n).SpikeRateOfChange = RateOfChange(PeakLocation(n)-PreSpikeTime : PeakLocation(n+1) + PostSpikeTime);
               SpikesExtracted=SpikesExtracted+1;
           else
               PW = sprintf('Warning the first spike was not extracted');
               disp(PW)
               SpikesNotExtracted = 1;
           end
       elseif n <NumSpikes 
           SpikesData.Spike(n).SpikeTime = Time(PeakLocation(n)-PreSpikeTime : PeakLocation(n+1) - PreSpikeTime);
           SpikesData.Spike(n).SpikeVoltage = Voltage(PeakLocation(n)-PreSpikeTime : PeakLocation(n+1) - PreSpikeTime);
           SpikesData.Spike(n).SpikeRateOfChange = RateOfChange(PeakLocation(n)-PreSpikeTime : PeakLocation(n+1) - PreSpikeTime);
           SpikesExtracted=SpikesExtracted+1;
       %This last spike can be extracted using the rest of Data
       %PROBLEM if the last spike ends early in the sweep the "nth spike"
       %will be very long since it ends with the sweep.
       elseif n == NumSpikes
           SpikesData.Spike(n).SpikeTime = Time(PeakLocation(n)-PreSpikeTime : end);
           SpikesData.Spike(n).SpikeVoltage = Voltage(PeakLocation(n)-PreSpikeTime : end);
           SpikesData.Spike(n).SpikeRateOfChange = RateOfChange(PeakLocation(n)-PreSpikeTime: end);
           SpikesExtracted=SpikesExtracted+1;
       end
   end
   
    P1=sprintf('Extracted %d Spikes', SpikesExtracted);
    disp(P1)
   if SpikesNotExtracted >= 1
        P2 =sprintf('Run Spike Deletion to delete spikes that were not extracted');
        disp(P2)
   end