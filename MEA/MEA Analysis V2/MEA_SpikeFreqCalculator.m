function [InstFreq , GloFreq] = MEA_SpikeFreqCalculator(PeakData,ChannelNumber,Time)
%This is meant to determine the spike frequency data
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 18, 2018
%LAST MODIFIED: April 19, 2018
%v1.0


%Check if there is Spike Data to load
if isempty(PeakData(ChannelNumber).SpikePeakTime)
    InstFreq = 0;
    GloFreq  = 0;
else
    %Load relevant data.
    SpikeTime = PeakData(ChannelNumber).SpikePeakTime;

    %take into consideration channels that lacked more than 1 spike, making
    %ISI not possible to calculate

    if length(SpikeTime) <= 1
    ISI = 0;
    elseif length(SpikeTime) > 1
        
        %Calculate the ISIs
        for ii = 1:(length(SpikeTime)-1)
            ISI(ii) =  SpikeTime(ii+1) - SpikeTime(ii);
        end
    end
    InstFreq = 1/mean(ISI);
    GloFreq  = PeakData(ChannelNumber).SpikeNum/Time;
end


