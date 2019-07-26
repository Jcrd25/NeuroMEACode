function MEA_SpectogramPlotter(MEAData, ChannelNumber, InterFileInterval,Stage)
%This function is meant to be used to calculate the power spectrum of the
%data adquired with the MEA multichannel systems. 
%
%Inputs
%   MEAData           - Data
%   InterFileInterval - Time in between files in ms
%   ChannelNumber     - Desired Channel Number to get the power spectrum
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 13, 2017
%LAST MODIFIED: November 3, 2017
%v1.0

SamplingRate = MEAData(1).RawData.SamplingRate;
% Check the setting to determine which MEA_FFT function is going to be used
for ii = 1:length(MEAData)
    TimePerRecording.File(ii) = length(MEAData(ii).RawData.VoltageTrace);
end
%Define the time variable 
TimeResolution = 10;%in seconds
TotalTime = 0;
for ii = 1:length(MEAData)
    TotalTime = TotalTime + TimePerRecording.File(ii);
end
TotalTime = (TotalTime + length(MEAData)*InterFileInterval)/SamplingRate;
%define the Time variable
Time = 1:TimeResolution:TotalTime;

%Calculate the power for all the files

n=1;
for ii = 1:length(MEAData)
    LengthData = length(MEAData(ii).RawData.VoltageTrace);
    %Calculate the Fast Fourier Transform
    %Desired window of analysis
    DesiredLength = TimeResolution*SamplingRate;
    %determine in how many segments you can dived the recording in order to
    %calculate the power of the time segment desired
    Fraction = floor(LengthData / (DesiredLength));
    Index1(1) = 1;
    Index2(1) = DesiredLength;
    for iii = 2:Fraction
        Index1(iii) = Index1(iii-1) + DesiredLength;
        Index2(iii) = Index2(iii-1) + DesiredLength;
    end
    %Calculate the Power for every time segment
    for iii = 1:Fraction
        VoltageTraceSeg = MEAData(ii).RawData.VoltageTrace(Index1(iii):Index2(iii),ChannelNumber);
        [Power(:,n), Frequency] = M_AverageFFT(VoltageTraceSeg, ChannelNumber, MEAData(ii).Settings.FileName, SamplingRate, 'N');
        n=n+1;
    end
end


%Plot the Spectogram
figure
imagesc((Time)*(1/60),Frequency,Power)
ylabel('Frequency (Hz)')
xlabel('Time (min)')
c = colorbar;
ylabel(c, 'dB/Hz')
colormap hot
ylim([1 100])
caxis([0 20]);
ham = gca;
set(ham, 'YDir', 'normal')
ChanNum = MEAData(1).Settings.InfoChannel.Label(ChannelNumber);
ChannelNumberDouble = str2double(ChanNum);
P1 = sprintf('Power Spectrum of Channel %d during %s', ChannelNumberDouble,Stage); 
title(P1)




