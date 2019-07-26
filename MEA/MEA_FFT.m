function MEA_FFT(VoltageTrace, ChannelNumber, FileName, SamplingRate)
%This funciton is meant to be used in order to obtain the periodrogram
%using the Fast Fourier Transform.
%
%Inputs
%   ChannelData - The voltage trace of the desired channel in µV
%   ChannelNumber - The number of the channel you want to obtain the
%   periodrogram
%   FileName - Name of the file. It is currently used in the tittle of the
%   periodogram
%   FrequencyUsed - Sampling rate used in the data adquisition
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
% FilteredChannelData = ChannelData;
%Substract the mean of the signal in order to get rid of the "0 Hz noise"
%if the mean is 0 then this function doesn't subtract anything
DFC_Data = VoltageTrace - mean(VoltageTrace);

%determine the length of the data
LengthData = length(DFC_Data);
%Calculate the Fast Fourier Transform
DFC_Data_dft = fft(DFC_Data);
%Take the positive values only
DFC_Data_dft = DFC_Data_dft(1:LengthData/2+1);
%Caculate the power
psd_DFC_Data = (1/(SamplingRate*LengthData)) * abs(DFC_Data_dft).^2;
psd_DFC_Data(2:end-1) = 2*psd_DFC_Data(2:end-1);
%define the frequency for x-axis
freq = 0:SamplingRate/length(DFC_Data):SamplingRate/2;
%Plot
figure
plot(freq,10*log10(psd_DFC_Data))
grid on
ChannelNumberDouble = str2double(ChannelNumber);
TitleName = M_Namer(FileName);
P1 = sprintf('Periodogram of Channel %d from file %s', ChannelNumberDouble, TitleName); 
title(P1)
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
%set the limits to better observed desired frequencies
xlim([0 100])
yl = ylim;
ylim([-50,yl(2)])
% figure
% plot(0:1/FrequencyUsed:((LengthofData-1)/FrequencyUsed), FilteredChannelData)
