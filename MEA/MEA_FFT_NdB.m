function MEA_FFT_NdB(ChannelData, ChannelNumber, FileName, SamplingRate)
%This funciton is meant to be used in order to obtain the periodrogram
%using the Fast Fourier Transform.
%
%Inputs
%   ChannelData - The voltage trace of the desired channel in µV
%   ChannelNumber - The number of the channel you want to obtain the
%   periodrogram
%   FileName - Name of the file. It is currently used in the tittle of the
%   periodogram
%   SamplingRate - Sampling rate used in the data adquisition
%Output
%   A periodrogram of the specified channel
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 14, 2017
%LAST MODIFIED: June 21, 2017
%v1.0
%NOTE: Modified from MEA_FFT and ForKevin


%Filter the data using a low pass Filter
%Filter Building tool

% %High Pass
% Fstop = 0.1;
% Fpass = 1;
% Astop = 65;
% Apass = 0.5;
% Fs    = SamplingRate;
% 
% 
% d = designfilt('highpassiir','StopbandFrequency',Fstop ,...
%   'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
%   'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','butter');
% 
% FilteredChannelData = filtfilt(d,ChannelData);

% %Low Pass
% Fpass = 400;
% Fstop = 500;
% Apass = 0.5;
% Astop = 65;
% Fs = SamplingRate;
% 
% d = designfilt('lowpassiir', ...
%   'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
%   'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
%   'DesignMethod','butter','SampleRate',Fs);

% %Bandpass filter
% Forder = 20;
% Freq1  = 20;
% Freq2  = 120;
% Fs     = SamplingRate;
% 
% d = designfilt('bandpassiir','FilterOrder',Forder, ...
%          'HalfPowerFrequency1',Freq1,'HalfPowerFrequency2',Freq2, ...
%          'SampleRate',Fs);

%Filter the data
% FilteredChannelData = filtfilt(d,ChannelData);
FilteredChannelData = ChannelData;
% FilteredChannelData = ChannelData;
%Substract the mean of the signal in order to get rid of the "0 Hz noise"
%if the mean is 0 then this function doesn't subtract anything
Demeaned_FilteredlData = FilteredChannelData - mean(FilteredChannelData);

%Define the sampling interval
dt = 1/SamplingRate;            

%Determine the time of the recording in s
LengthofData = length(Demeaned_FilteredlData);
RecordingTime = LengthofData/SamplingRate;

%Calculate the Fast Fourier Transform

Sxx = 2*dt^2/RecordingTime * fft(Demeaned_FilteredlData).*conj(fft(Demeaned_FilteredlData));              % Compute the power spectrum. 
Sxx = Sxx(1:length(Demeaned_FilteredlData)/2+1);                         % ... ignore negative frequencies.	
df  = 1/RecordingTime;                                           % Determine the frequency resolution. 
fNQ = 1/dt/2;                                         % Determine the Nyquist frequency. 
faxis = (0:df:fNQ);                                 % Construct the frequency axis.
figure
plot(faxis, Sxx)                                    % Plot power versus frequency. 
yl = ylim;
ylim([0 ,yl(2)])                                       % Set range of y-values in plot.
xlim([0 100])
xlabel('Frequency [Hz]'); ylabel('Power [\muV^2/Hz]')
ChannelNumberDouble = str2double(ChannelNumber);
TitleName = M_Namer(FileName);
P1 = sprintf('Power Spectrum of Channel %d from file %s', ChannelNumberDouble, TitleName); 
title(P1)

