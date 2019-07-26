function MEA_FFT_Spectrogram(ChannelData, ChannelNumber, FileName, SamplingRate)
%This funciton is meant to be used in order to obtain the spectrogram
%using the Fast Fourier Transform.
%
%Inputs
%   ChannelData   - The voltage trace of the desired channel in µV
%   ChannelNumber - The number of the channel you want to obtain the
%                   periodrogram
%   FileName      - Name of the file. It is currently used in the tittle of the
%                   periodogram
%   SamplingRate  - Sampling rate used in the data adquisition
%Output
%   A periodrogram of the specified channel
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 14, 2017
%LAST MODIFIED: June 12, 2017
%v1.0



% % %Filter the data using a low pass Filter
% %Design the low pass filter
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
% 

% %Bandpass filter
% Forder = 20;
% Freq1  = 1;
% Freq2  = 120;
% Fs     = SamplingRate;
% 
% d = designfilt('bandpassiir','FilterOrder',Forder, ...
%          'HalfPowerFrequency1',Freq1,'HalfPowerFrequency2',Freq2, ...
%          'SampleRate',Fs);
%      
%Filter the data
% FilteredChannelData = filtfilt(d,ChannelData);
FilteredChannelData = ChannelData;
%Substract the mean of the signal in order to get rid of the "0 Hz noise"
%if the mean is 0 then this function doesn't subtract anything
Demeaned_FilteredlData = FilteredChannelData - mean(FilteredChannelData);

%Use a Hamming Window that advances 100 ms and calculate the FFT in 
%bins. Seperate the time in bins defeined by their middle value

%For gamma oscillations Dr. Omar Ahmed recommended to use  half a second
%bins and an advance of 100 ms.

BinsSec     = 2 ; %in s
DesiredRoll = 100; %in ms

%Convert DesiredRoll to s
Rolling = DesiredRoll*0.001;
%Calculate the Total Time in the recording
TotalTime = floor(length(Demeaned_FilteredlData)/SamplingRate);
%Calculate the middle time value for each bin
tmid = (0+(BinsSec/2)):Rolling:(TotalTime-(BinsSec/2));
%Define the frequency                      
Fs = SamplingRate;
spec = [];
%Define the indexes used to generate the bins
Index = 1:(Rolling*SamplingRate):(length(Demeaned_FilteredlData));

for ii = 1:(length(tmid)-1)
    %Select the data segment that is going to be used to calculate the FFT
    LeftEnd  =  Index(ii);
    RightEnd =  Index(ii+(BinsSec/Rolling));
    x = Demeaned_FilteredlData(LeftEnd:RightEnd);
    %Calculate the FFT for the desired time segment
    
    %Calculate the hamming window
    N = length(x);
    ham = hamming(N);
    %Multiply the data by the hamming window
    hx = ham.*x;
    %Perform the FFT on the data
    xdft = fft(hx);
    xdft = xdft(1:round(N/2)+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(hx):Fs/2;
    spec(:,ii) = 10*log10(psdx);
end

%Plot the Spectogram using the decibel scale
figure
imagesc(tmid,freq,spec)
ylabel('Frequency (Hz)')
xlabel('Time (s)')
c = colorbar;
ylabel(c, 'dB/Hz')
colormap hot
ylim([1 100])
caxis([-10 20]);
ham = gca;
set(ham, 'YDir', 'normal')
ChannelNumberDouble = str2double(ChannelNumber);
TitleName = M_Namer(FileName);
P1 = sprintf('Power Spectrum of Channel %d from file %s', ChannelNumberDouble, TitleName); 
title(P1)

% %This is a spectogram using the wavelet transform
% figure
% cwt(Demeaned_FilteredlData , SamplingRate)

% figure
%     %Plot the spectogram
%     spectrogram(Demeaned_FilteredlData, 256*100 , [], [] , SamplingRate,'yaxis');
%      ylim([0 0.1])
%      caxis([-8 7]);
%      colormap hot
%      %Put title to the graph
% %      ChannelNumberDouble = str2double(ChannelNumber);
% %      P1 = sprintf('Spetogram of Channel %d from file %s', ChannelNumberDouble, FileName); 
% %      title(P1)

