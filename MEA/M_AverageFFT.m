function [Power, Frequency] = M_AverageFFT(VoltageTrace, ChannelNumber, FileName, SamplingRate, PlotAnswer,Color)
%This funciton is meant to be used in order to obtain the periodrogram
%using the Fast Fourier Transform by averaging the trace
%
%Inputs
%   VoltageTrace  - The voltage trace of the desired channel in µV
%   ChannelNumber - The number of the channel you want to obtain the
%   periodrogram
%   FileName      - Name of the file. It is currently used in the tittle of the
%   periodogram
%   FrequencyUsed - Sampling rate used in the data adquisition
%   PlotAnswer    - What setting to use to plot the figures
%       N - Does not plot
%       Y - Plots every one individually
%       S - Plots them in one figure
%Output
%   A periodrogram of the specified channel
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 12, 2017
%LAST MODIFIED: June 21, 2017
%v1.0

%Filter the data using a low pass Filter
% Low pass filter
%Design the low pass filter

Fpass = 400;
Fstop = 500;
Apass = 0.5;
Astop = 65;
Fs = SamplingRate;

d = designfilt('lowpassiir', ...
  'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
  'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
  'DesignMethod','butter','SampleRate',Fs);

VoltageTrace = filtfilt(d,VoltageTrace);

%Substract the mean of the signal in order to get rid of the "0 Hz noise"
%if the mean is 0 then this function doesn't subtract anything
VoltageTrace = VoltageTrace - mean(VoltageTrace);

%determine the length of the data
LengthData = length(VoltageTrace);
%Calculate the Fast Fourier Transform
%Desired window of analysis
DesiredTimeLength = 1; %in s
DesiredLength = DesiredTimeLength*SamplingRate;
FFT_Data = struct;
%determine in how many segments you can dived the recording in order to
%calculate the average
Fraction = floor(LengthData / (DesiredLength));

Index1(1) = 1;
Index2(1) = DesiredLength;
for ii = 2:Fraction
    Index1(ii) = Index1(ii-1) + DesiredLength;
    Index2(ii) = Index2(ii-1) + DesiredLength;
end
for ii = 1:Fraction
    [FFT_Data(ii).Power, Frequency] = M_FFT(VoltageTrace(Index1(ii):Index2(ii)), SamplingRate);
end
Power = zeros(Fraction,length(FFT_Data(1).Power));
for n = 1:Fraction
    for ii = 1:length(FFT_Data(1).Power)
    Power(n,ii) = (FFT_Data(n).Power(ii));
    end
end
Power = mean(Power);    

if ~exist('PlotAnswer','var')
     %if PlotAnswer Setting does not exist, so default it all
      PlotAnswer = 'Y';
      fprintf('No Setting specified.\n Using the default setting\n')
end

if PlotAnswer == 'N'
elseif PlotAnswer == 'S'
    %Plot
%     figure
%     hold on
    plot(Frequency, 10*log10(Power),Color)
    grid on
    ChannelNumberDouble = str2double(ChannelNumber);
%     TitleName = M_Namer(FileName);
%     P1 = sprintf('Periodogram of Channel %d from file %s', ChannelNumberDouble, TitleName); 
    P1 = sprintf('Periodogram of Channel %d', ChannelNumberDouble); 
    title(P1)
    xlabel('Frequency (Hz)')
    ylabel('Power/Frequency (dB/Hz)')
    %set the limits to better observed desired frequencies
    xlim([0 100])
    
elseif PlotAnswer == 'Y'

    %Plot
    figure
    plot(Frequency,10*log10(Power),Color)
    grid on
    ChannelNumberDouble = str2double(ChannelNumber);
    TitleName = M_Namer(FileName);
    P1 = sprintf('Periodogram of Channel %d from file %s', ChannelNumberDouble, TitleName); 
    title(P1)
    xlabel('Frequency (Hz)')
    ylabel('Power/Frequency (dB/Hz)')
    %set the limits to better observed desired frequencies
    xlim([0 100])
end

