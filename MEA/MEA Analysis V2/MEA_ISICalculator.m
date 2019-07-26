function ISIData = MEA_ISICalculator(PeakData,ChannelNumber)
%This is meant to determine the interspike interval.
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 12, 2018
%LAST MODIFIED: April 18, 2018
%v1.0


%Check if there is Spike Data to load
if isempty(PeakData(ChannelNumber).SpikePeakTime)
    ISI = 0;
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
end
%Save Data
ISIData = ISI';